"""Pointwise kinetics for Darcy + diluted-species DFM (model sheet, SI units)."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

from dfm_config import DfmConfig

from fenicsx.properties import (
    R_gas,
    T_isothermal_K,
    molar_mass_kg,
    partial_pressure,
    q_e_const_mol_kg,
    rho_a_bed_kg_m3,
    rho_a_over_epsilon,
    species_list,
)


@dataclass(frozen=True)
class KineticRates:
    """Solid-phase rates [mol/(kg_s·s)]: r_ads, r_des, r_CH4."""

    r_ads: float = 0.0
    r_des: float = 0.0
    r_ch4: float = 0.0


def _csc_params(cfg: DfmConfig) -> tuple[float, float, float, float, float]:
    """Return (q_e_scale, K_C, K_R, K_c, a). q_e_scale is mandatory COMSOL CSC q_e."""
    ads = cfg.solid_kinetics.get("adsorption_CO2") or {}
    csc = ads.get("csc") or {}
    q_e = csc.get("q_e_mol_kg")
    if q_e is None or float(q_e) <= 0:
        raise ValueError("adsorption_CO2.csc.q_e_mol_kg > 0 is required (CSC capacity scale).")
    K_C = float(csc.get("K_C_Pa_inv") or 0.0)
    K_R = float(csc.get("K_R", 0.0))
    K_c = float(csc.get("K_c_Pa_inv", K_C))
    a = float(csc.get("a", 1.0))
    if K_C <= 0:
        raise ValueError("adsorption_CO2.csc.K_C_Pa_inv required.")
    return float(q_e), K_C, K_R, K_c, a


_ATM_PA = 101325.0


def _csc_pressure_scale(cfg: DfmConfig) -> float:
    """Pa → CSC argument unit (COMSOL CSC declares argunit=atm; K_C is per that unit)."""
    ads = cfg.solid_kinetics.get("adsorption_CO2") or {}
    csc = ads.get("csc") or {}
    unit = str(csc.get("P_unit", "Pa")).strip().lower()
    if unit == "atm":
        return 1.0 / _ATM_PA
    if unit == "bar":
        return 1e-5
    return 1.0


def q_e_csc(P_CO2_Pa: float, cfg: DfmConfig) -> float:
    """
    CSC isotherm q_e(P_CO2):

        q_e(P) = q_e_scale * K_C P (1 + (a+1) K_R P^a) / (1 + K_C P + K_c K_R P^(a+1))

    with q_e_scale = csc.q_e_mol_kg (mandatory). P is converted from Pa to csc.P_unit
    (COMSOL POW_V07: atm) before evaluation.
    """
    q_e_scale, K_C, K_R, K_c, a = _csc_params(cfg)
    P = max(float(P_CO2_Pa), 0.0) * _csc_pressure_scale(cfg)
    Pa = P**a
    Pap1 = P ** (a + 1.0)
    num = q_e_scale * K_C * P * (1.0 + (a + 1.0) * K_R * Pa)
    den = 1.0 + K_C * P + K_c * K_R * Pap1
    return max(num / max(den, 1e-30), 0.0)


def q_e_equilibrium(C_CO2: float, cfg: DfmConfig) -> float:
    return q_e_csc(partial_pressure(C_CO2, cfg), cfg)


def theta_from_q(q: float, cfg: DfmConfig) -> float:
    return max(float(q), 0.0) / q_e_const_mol_kg(cfg)


def _k_ads_si(ads: dict[str, Any]) -> float:
    k = ads.get("k_ads")
    if k is None:
        raise ValueError("adsorption_CO2.k_ads required.")
    k = float(k)
    unit = (ads.get("k_ads_unit") or "1/s").strip()
    if unit in ("1/min", "min^-1", "per_min"):
        return k / 60.0
    return k


def _k_des_si(cfg: DfmConfig) -> float:
    """
    Desorption coefficient for

        r_des = k_des · C_H2 · θ   [mol/(kg_solid·s)]

    Dimensional consistency (POW / COMSOL):
        [C_H2] = mol/m³,  [θ] = 1  ⇒  [k_des] = m³/(kg·s).

    Accepted units (converted to m³/(kg·s)):
        m3/(kg*s) and aliases — used as-is;
        cm3/(g*min) — × 1/60000.
    ``1/s`` is rejected.
    """
    ads = cfg.solid_kinetics.get("adsorption_CO2") or {}
    k = ads.get("k_des")
    if k is None:
        return 0.0
    k = float(k)
    unit = str(ads.get("k_des_unit") or "m3/(kg*s)").strip().lower().replace(" ", "")
    # Normalize unicode superscripts / middots for matching.
    unit = (
        unit.replace("³", "3")
        .replace("^3", "3")
        .replace("·", "*")
        .replace("⋅", "*")
    )
    si_aliases = {
        "m3/(kg*s)",
        "m3_kg_s",
        "m3/kg/s",
        "m3/(kg.s)",
    }
    cm3_g_min = {
        "cm3/(g*min)",
        "cm3/(g.min)",
        "cm3/g/min",
        "cm3_g_min",
    }
    if unit in si_aliases or unit == "":
        return k
    if unit in cm3_g_min:
        return k / 60000.0
    if unit in ("1/s", "s^-1", "s-1", "1/min", "min^-1", "per_min"):
        raise ValueError(
            f"adsorption_CO2.k_des_unit={ads.get('k_des_unit')!r} is invalid for "
            "r_des = k_des·C_H2·θ. Use m3/(kg*s) or cm3/(g*min)."
        )
    raise ValueError(
        f"Unsupported adsorption_CO2.k_des_unit={ads.get('k_des_unit')!r}; "
        "expected m3/(kg*s) (or cm3/(g*min))."
    )


def _meth_loading_driver(cfg: DfmConfig) -> str:
    meth = cfg.solid_kinetics.get("methanation") or {}
    arr = meth.get("arrhenius") or {}
    driver = arr.get("loading_driver") or meth.get("loading_driver") or "theta"
    return str(driver).strip()


def _meth_loading_value(q: float, cfg: DfmConfig) -> float:
    qv = max(float(q), 0.0)
    if _meth_loading_driver(cfg) in ("q_CO2", "q", "qCO2"):
        return qv
    return theta_from_q(qv, cfg)


def _norm_unit(unit: str | None) -> str:
    return (
        str(unit or "")
        .strip()
        .lower()
        .replace(" ", "")
        .replace("·", "*")
        .replace("⋅", "*")
        .replace("³", "3")
        .replace("^3", "3")
    )


def _a_rxn_si(arr: dict[str, Any]) -> float:
    """
    Pre-exponential A_rxn → mol/(kg·s) (same convention as r_CH4; m1/m2 ignored in the unit).

    Accepted: mol/(kg*s) (default), kmol/(kg*s) → ×1000.
    """
    A = arr.get("A_rxn")
    if A is None:
        raise ValueError("methanation.arrhenius.A_rxn required.")
    A = float(A)
    unit = _norm_unit(arr.get("A_rxn_unit") or "mol/(kg*s)")
    mol_aliases = {"mol/(kg*s)", "mol/kg/s", "mol/(kg.s)", "mol_kg_s", ""}
    kmol_aliases = {"kmol/(kg*s)", "kmol/kg/s", "kmol/(kg.s)", "kmol_kg_s"}
    if unit in mol_aliases:
        return A
    if unit in kmol_aliases:
        return A * 1000.0
    raise ValueError(
        f"Unsupported methanation.arrhenius.A_rxn_unit={arr.get('A_rxn_unit')!r}; "
        "expected mol/(kg*s) or kmol/(kg*s)."
    )


def _ea_si_J_mol(arr: dict[str, Any]) -> float:
    """Activation energy → J/mol. Accepts Ea_J_mol (+ optional Ea_unit)."""
    if "Ea_J_mol" in arr and arr.get("Ea_J_mol") is not None:
        Ea = float(arr["Ea_J_mol"])
    elif arr.get("Ea") is not None:
        Ea = float(arr["Ea"])
    else:
        return 0.0
    unit = _norm_unit(arr.get("Ea_unit") or "J/mol")
    if unit in ("j/mol", "j·mol^-1", "j_mol", ""):
        return Ea
    if unit in ("kj/mol", "kj·mol^-1", "kj_mol"):
        return Ea * 1000.0
    raise ValueError(
        f"Unsupported methanation.arrhenius.Ea_unit={arr.get('Ea_unit')!r}; "
        "expected J/mol or kJ/mol."
    )


def _meth_exponents(cfg: DfmConfig) -> tuple[float, float, float]:
    meth = cfg.solid_kinetics.get("methanation") or {}
    arr = meth.get("arrhenius") or {}
    m1 = float(arr.get("m1", 1.0))
    m2 = float(arr.get("m2", 1.0))
    A = _a_rxn_si(arr)
    Ea = _ea_si_J_mol(arr)
    T = T_isothermal_K(cfg)
    k_rxn = float(A)
    if Ea != 0.0:
        k_rxn *= math.exp(-Ea / (R_gas(cfg) * T))
    return k_rxn, m1, m2


def _meth_ph2_scale(cfg: DfmConfig) -> float:
    """Pa → PH2 unit of the rate law (COMSOL kinetics() declares PH2 in bar)."""
    meth = cfg.solid_kinetics.get("methanation") or {}
    arr = meth.get("arrhenius") or {}
    unit = str(arr.get("PH2_unit", "Pa")).strip().lower()
    if unit == "bar":
        return 1e-5
    if unit == "atm":
        return 1.0 / _ATM_PA
    if unit == "kpa":
        return 1e-3
    return 1.0


def kinetic_rates(
    *,
    C: dict[str, float],
    q: float,
    cfg: DfmConfig,
    s_co2: float,
    s_h2: float,
) -> KineticRates:
    """
    Model-sheet kinetics (gates s_CO2, s_H2):

        r_ads = k_ads (q_e − q) s_CO2
        r_CH4 = k_rxn P_H2^m1 · (q or θ)^m2 s_H2   (loading_driver: q_CO2 | theta)
        r_des = k_des C_H2 θ s_H2  (≡ r_CO2 in COMSOL POW_V07)
    """
    ads_cfg = cfg.solid_kinetics.get("adsorption_CO2") or {}
    qv = max(float(q), 0.0)
    theta = theta_from_q(qv, cfg)
    load = _meth_loading_value(q, cfg)
    c_h2 = max(float(C.get("H2", 0.0)), 0.0)
    p_h2 = partial_pressure(c_h2, cfg) * _meth_ph2_scale(cfg)

    r_ads = 0.0
    if s_co2 > 0:
        qe = q_e_equilibrium(C["CO2"], cfg)
        r_ads = _k_ads_si(ads_cfg) * (qe - qv) * s_co2

    r_des = 0.0
    r_ch4 = 0.0
    if s_h2 > 0:
        k_des = _k_des_si(cfg)
        r_des = k_des * c_h2 * max(theta, 0.0) * s_h2
        k_rxn, m1, m2 = _meth_exponents(cfg)
        r_ch4 = k_rxn * (max(p_h2, 0.0) ** m1) * (max(load, 0.0) ** m2) * s_h2

    return KineticRates(r_ads=r_ads, r_des=r_des, r_ch4=r_ch4)


def _stoichiometric_molar_terms(kr: KineticRates) -> dict[str, float]:
    """
    Stoichiometric factors on solid rates r [mol/(kg_s·s)] — before volume prefactor.

    Shared by Darcy (× ρ_a M_i) and diluted transport (× ρ_a/ε).
    """
    return {
        "CO2": -kr.r_ads + kr.r_des,
        "H2": -4.0 * kr.r_ch4,
        "CH4": kr.r_ch4,
        "H2O": 2.0 * kr.r_ch4,
        "N2": 0.0,
    }


def dq_dt_solid(
    *,
    C: dict[str, float],
    q: float,
    cfg: DfmConfig,
    s_co2: float,
    s_h2: float,
) -> float:
    """dq/dt = r_ads − r_CH4 − r_des  (equivalent to ∂θ/∂t = (r_ads − r_CH4 − r_des)/q_{e,const})."""
    kr = kinetic_rates(C=C, q=q, cfg=cfg, s_co2=s_co2, s_h2=s_h2)
    return kr.r_ads - kr.r_ch4 - kr.r_des


def mass_source_sum_Ri(
    *,
    C: dict[str, float],
    q: float,
    cfg: DfmConfig,
    s_co2: float,
    s_h2: float,
) -> float:
    """
    Σ R_i [kg/(m³_bed·s)] for Darcy continuity (slide: mass basis, bed volume).

        R_i = ρ_a M_i × (stoichiometric term on r)
        Σ R_i = ρ_a [ M_CO2(-r_ads+r_des) + M_H2(-4 r_CH4) + M_CH4 r_CH4 + M_H2O(2 r_CH4) ]

    ρ_a = (1−ε) ρ_s — no factor 1/ε; molar masses convert mol→kg.
    """
    rho_a = rho_a_bed_kg_m3(cfg)
    kr = kinetic_rates(C=C, q=q, cfg=cfg, s_co2=s_co2, s_h2=s_h2)
    terms = _stoichiometric_molar_terms(kr)
    return rho_a * sum(molar_mass_kg(sp, cfg) * terms[sp] for sp in terms)


def k_ads_si(cfg: DfmConfig) -> float:
    """Adsorption LDF constant k_ads [1/s] (unit-converted from config)."""
    return _k_ads_si(cfg.solid_kinetics.get("adsorption_CO2") or {})


def csc_qe_evaluator(cfg: DfmConfig):
    """
    Fast closure q_e(C_CO2 [mol/m³]) → [mol/kg] with precomputed CSC constants.

    Accepts scalars or numpy arrays (vectorized).
    """
    q_e_scale, K_C, K_R, K_c, a = _csc_params(cfg)
    RT_scale = R_gas(cfg) * T_isothermal_K(cfg) * _csc_pressure_scale(cfg)

    def qe(C_mol_m3):
        P = np.maximum(C_mol_m3, 0.0) * RT_scale
        Pa = P**a
        num = q_e_scale * K_C * P * (1.0 + (a + 1.0) * K_R * Pa)
        den = 1.0 + K_C * P + K_c * K_R * Pa * P
        return np.maximum(num / np.maximum(den, 1e-30), 0.0)

    return qe


def adsorption_exchange_point(
    C0: float,
    q0: float,
    dt: float,
    *,
    k_ads: float,
    pref: float,
    qe,
) -> tuple[float, float]:
    """
    Backward-Euler step of the local LDF exchange (mass-conserving):

        dq/dt = k_ads (q_e(C) − q),   dC/dt = −(ρ_a/ε) dq/dt
        ⇒ C + (ρ_a/ε) q = const

    Reduces to a monotone scalar equation in q_new; solved by bisection.
    Replaces the explicit r_ads transport source, which is unstable/lossy for
    dt ≫ ε/(ρ_a k_ads dq_e/dC).
    """
    C0 = max(float(C0), 0.0)
    q0 = max(float(q0), 0.0)
    if dt <= 0.0 or k_ads <= 0.0:
        return C0, q0
    qe0 = qe(C0)
    if qe0 > q0:  # adsorption: q ∈ [q0, q0 + C0/pref]
        lo, hi = q0, q0 + C0 / pref
    else:  # LDF desorption: q ∈ [0, q0]
        lo, hi = 0.0, q0
    if hi - lo < 1e-15:
        return C0, q0

    def g(qn: float) -> float:
        Cn = C0 - pref * (qn - q0)
        return qn - q0 - dt * k_ads * (qe(max(Cn, 0.0)) - qn)

    if g(lo) >= 0.0:
        return C0, q0
    if g(hi) <= 0.0:
        qn = hi
    else:
        for _ in range(52):
            mid = 0.5 * (lo + hi)
            if g(mid) <= 0.0:
                lo = mid
            else:
                hi = mid
        qn = 0.5 * (lo + hi)
    Cn = max(C0 - pref * (qn - q0), 0.0)
    return Cn, qn


def make_reaction_exchange(cfg: DfmConfig):
    """
    Pointwise backward-Euler step for methanation + H2-driven desorption
    (mass-conserving Lie splitting; replaces explicit r_CH4/r_des sources):

        x ≡ dt·r_CH4:  H_new = H0 − 4 (ρ_a/ε) x
                       q_new (1 + dt k_des H_new/q_ref) = q0 − x   [closed form]
        residual F(x) = x − dt k_rxn P_H2(H_new)^m1 load(q_new)^m2 = 0  (monotone)

    Outputs (C_H2', q', dCH4, dH2O, dCO2_des) per node; bisection on x.
    """
    k_rxn, m1, m2 = _meth_exponents(cfg)
    k_des = _k_des_si(cfg)
    pref = rho_a_over_epsilon(cfg)
    q_ref = q_e_const_mol_kg(cfg)
    RT_ph2 = R_gas(cfg) * T_isothermal_K(cfg) * _meth_ph2_scale(cfg)
    driver_is_q = _meth_loading_driver(cfg) in ("q_CO2", "q", "qCO2")

    def step(H0: float, q0: float, dt: float) -> tuple[float, float, float, float, float]:
        H0 = max(float(H0), 0.0)
        q0 = max(float(q0), 0.0)
        if dt <= 0.0 or H0 <= 0.0 or q0 <= 0.0:
            return H0, q0, 0.0, 0.0, 0.0

        def q_of_x(x: float) -> float:
            Hn = max(H0 - 4.0 * pref * x, 0.0)
            return max((q0 - x) / (1.0 + dt * k_des * Hn / q_ref), 0.0), Hn

        def F(x: float) -> float:
            qn, Hn = q_of_x(x)
            if k_rxn <= 0.0:
                return x
            load = qn if driver_is_q else qn / q_ref
            r = k_rxn * (Hn * RT_ph2) ** m1 * load**m2
            return x - dt * r

        x_hi = min(H0 / (4.0 * pref), q0)
        if x_hi <= 0.0 or F(0.0) >= 0.0:
            x = 0.0
        elif F(x_hi) <= 0.0:
            x = x_hi
        else:
            lo, hi = 0.0, x_hi
            for _ in range(52):
                mid = 0.5 * (lo + hi)
                if F(mid) <= 0.0:
                    lo = mid
                else:
                    hi = mid
            x = 0.5 * (lo + hi)

        qn, Hn = q_of_x(x)
        x_des = max(q0 - x - qn, 0.0)  # dt·r_des consistent with the implicit solve
        return Hn, qn, pref * x, 2.0 * pref * x, pref * x_des

    return step


_N_BISECT = 40  # ~1e-12 relative bracket width; enough for double


def adsorption_exchange_arrays(
    C: np.ndarray,
    q: np.ndarray,
    dt: float,
    *,
    k_ads: float,
    pref: float,
    qe,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Vectorized adsorption_exchange_point over nodal arrays (same semantics).

    Only nodes that need a solve are bisected (indexed subset), so idle /
    equilibrium nodes cost O(1) rather than 40× qe evaluations.
    """
    C0 = np.maximum(np.asarray(C, dtype=float), 0.0)
    q0 = np.maximum(np.asarray(q, dtype=float), 0.0)
    if dt <= 0.0 or k_ads <= 0.0:
        return C0, q0
    qe0 = qe(C0)
    # Skip when driving force is negligible everywhere (dead time / plateau)
    if float(np.max(np.abs(qe0 - q0))) * dt * k_ads < 1e-14:
        return C0, q0
    ads = qe0 > q0
    lo = np.where(ads, q0, 0.0)
    hi = np.where(ads, q0 + C0 / pref, q0)
    span = hi - lo
    cand = span >= 1e-15
    if not np.any(cand):
        return C0, q0

    def g_idx(qn: np.ndarray, idx: np.ndarray) -> np.ndarray:
        Cn = np.maximum(C0[idx] - pref * (qn - q0[idx]), 0.0)
        return qn - q0[idx] - dt * k_ads * (qe(Cn) - qn)

    idx_c = np.flatnonzero(cand)
    g_lo = g_idx(lo[idx_c], idx_c)
    active_c = g_lo < 0.0
    if not np.any(active_c):
        return C0, q0
    idx = idx_c[active_c]
    g_hi = g_idx(hi[idx], idx)
    take_hi = g_hi <= 0.0
    bis = ~take_hi

    qn = q0.copy()
    qn[idx[take_hi]] = hi[idx[take_hi]]
    if np.any(bis):
        ib = idx[bis]
        lo_b = lo[ib].copy()
        hi_b = hi[ib].copy()
        for _ in range(_N_BISECT):
            mid = 0.5 * (lo_b + hi_b)
            le = g_idx(mid, ib) <= 0.0
            lo_b = np.where(le, mid, lo_b)
            hi_b = np.where(le, hi_b, mid)
        qn[ib] = 0.5 * (lo_b + hi_b)

    Cn = C0.copy()
    Cn[idx] = np.maximum(C0[idx] - pref * (qn[idx] - q0[idx]), 0.0)
    return Cn, qn


def make_reaction_exchange_arrays(cfg: DfmConfig):
    """Vectorized make_reaction_exchange over nodal arrays (same semantics)."""
    k_rxn, m1, m2 = _meth_exponents(cfg)
    k_des = _k_des_si(cfg)
    pref = rho_a_over_epsilon(cfg)
    q_ref = q_e_const_mol_kg(cfg)
    RT_ph2 = R_gas(cfg) * T_isothermal_K(cfg) * _meth_ph2_scale(cfg)
    driver_is_q = _meth_loading_driver(cfg) in ("q_CO2", "q", "qCO2")

    def step(
        H: np.ndarray, q: np.ndarray, dt: float
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        H0 = np.maximum(np.asarray(H, dtype=float), 0.0)
        q0 = np.maximum(np.asarray(q, dtype=float), 0.0)
        zeros = np.zeros_like(H0)
        if dt <= 0.0:
            return H0, q0, zeros, zeros, zeros
        pos = (H0 > 0.0) & (q0 > 0.0)
        if not np.any(pos):
            return H0, q0, zeros, zeros, zeros

        def q_of_x_idx(x: np.ndarray, idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            Hn = np.maximum(H0[idx] - 4.0 * pref * x, 0.0)
            qn = np.maximum((q0[idx] - x) / (1.0 + dt * k_des * Hn / q_ref), 0.0)
            return qn, Hn

        def F_idx(x: np.ndarray, idx: np.ndarray) -> np.ndarray:
            qn, Hn = q_of_x_idx(x, idx)
            if k_rxn <= 0.0:
                return x
            load = qn if driver_is_q else qn / q_ref
            r = k_rxn * (Hn * RT_ph2) ** m1 * load**m2
            return x - dt * r

        idx_p = np.flatnonzero(pos)
        x_hi_p = np.minimum(H0[idx_p] / (4.0 * pref), q0[idx_p])
        cand = x_hi_p > 0.0
        if not np.any(cand):
            return H0, q0, zeros, zeros, zeros
        idx = idx_p[cand]
        x_hi = x_hi_p[cand]
        react = F_idx(np.zeros_like(x_hi), idx) < 0.0
        if not np.any(react):
            return H0, q0, zeros, zeros, zeros
        idx_r = idx[react]
        x_hi_r = x_hi[react]
        f_hi = F_idx(x_hi_r, idx_r)
        take_hi = f_hi <= 0.0
        bis = ~take_hi

        x_full = zeros.copy()
        x_full[idx_r[take_hi]] = x_hi_r[take_hi]
        if np.any(bis):
            ib = idx_r[bis]
            lo_b = np.zeros(ib.size)
            hi_b = x_hi_r[bis].copy()
            for _ in range(_N_BISECT):
                mid = 0.5 * (lo_b + hi_b)
                le = F_idx(mid, ib) <= 0.0
                lo_b = np.where(le, mid, lo_b)
                hi_b = np.where(le, hi_b, mid)
            x_full[ib] = 0.5 * (lo_b + hi_b)

        idx_all = idx_r
        qn_r, Hn_r = q_of_x_idx(x_full[idx_all], idx_all)
        Hn = H0.copy()
        qn = q0.copy()
        Hn[idx_all] = Hn_r
        qn[idx_all] = qn_r
        x = x_full
        x_des = np.maximum(q0 - x - qn, 0.0)
        d_ch4 = pref * x
        return Hn, qn, d_ch4, 2.0 * d_ch4, pref * x_des

    return step


def make_mass_source_evaluator(cfg: DfmConfig):
    """
    Vectorized Σ R_i [kg/(m³_bed·s)] for Darcy continuity — same math as
    mass_source_sum_Ri evaluated on nodal arrays.
    """
    rho_a = rho_a_bed_kg_m3(cfg)
    k_ads = k_ads_si(cfg)
    k_des = _k_des_si(cfg)
    k_rxn, m1, m2 = _meth_exponents(cfg)
    RT_ph2 = R_gas(cfg) * T_isothermal_K(cfg) * _meth_ph2_scale(cfg)
    q_ref = q_e_const_mol_kg(cfg)
    driver_is_q = _meth_loading_driver(cfg) in ("q_CO2", "q", "qCO2")
    qe = csc_qe_evaluator(cfg)
    M = {sp: molar_mass_kg(sp, cfg) for sp in ("CO2", "H2", "CH4", "H2O")}

    def sum_Ri(
        C_co2: np.ndarray,
        C_h2: np.ndarray,
        q: np.ndarray,
        s_co2: float,
        s_h2: float,
    ) -> np.ndarray:
        qv = np.maximum(np.asarray(q, dtype=float), 0.0)
        out = np.zeros_like(qv)
        if s_co2 > 0:
            r_ads = k_ads * (qe(C_co2) - qv) * s_co2
            out += M["CO2"] * (-r_ads)
        if s_h2 > 0:
            theta = qv / q_ref
            c_h2 = np.maximum(np.asarray(C_h2, dtype=float), 0.0)
            r_des = k_des * c_h2 * theta * s_h2
            load = qv if driver_is_q else theta
            r_ch4 = k_rxn * (c_h2 * RT_ph2) ** m1 * load**m2 * s_h2
            out += (
                M["CO2"] * r_des
                + M["H2"] * (-4.0 * r_ch4)
                + M["CH4"] * r_ch4
                + M["H2O"] * (2.0 * r_ch4)
            )
        return rho_a * out

    return sum_Ri


def gas_molar_sources_transport(
    *,
    C: dict[str, float],
    q: float,
    cfg: DfmConfig,
    s_co2: float,
    s_h2: float,
) -> dict[str, float]:
    """
    R_i [mol/(m³_void·s)] for diluted-species transport (slide: void volume).

        R_CO2  = (ρ_a/ε)(-r_ads + r_des)
        R_H2   = (ρ_a/ε)(-4 r_CH4)
        R_CH4  = (ρ_a/ε)(r_CH4)
        R_H2O  = (ρ_a/ε)(2 r_CH4)
        R_N2   = 0

    Prefactor ρ_a/ε — no molar mass M_i (unlike Darcy Σ R_i).
    """
    pref = rho_a_over_epsilon(cfg)
    kr = kinetic_rates(C=C, q=q, cfg=cfg, s_co2=s_co2, s_h2=s_h2)
    terms = _stoichiometric_molar_terms(kr)
    return {sp: pref * terms[sp] for sp in species_list(cfg)}
