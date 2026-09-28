"""
Wilke mixture diffusivity + COMSOL DL axial dispersion (v1.1, no COMSOL runtime).

COMSOL comp2 (var15):
  Dm_i = (1-w_i) / sum_{j!=i} ( y_j / Diff_{ij} )
  D_z  = 0.7 * Dm_i + 0.5 * dp * Vel_D   (Vel_D superficial)

Binary Diff_ij:
  fuller — Fuller et al. correlation at T_isothermal, P_total (runtime)
  user   — transport.binary_diffusion_m2_s table (values you supply)
"""

from __future__ import annotations

import math

import numpy as np

from dfm_config import DfmConfig
from fenicsx.properties import T_isothermal_K, molar_mass_kg

# COMSOL species index order in Diff_ij keys (1=N2 … 5=H2O).
WILKE_ORDER = ["N2", "CO2", "H2", "CH4", "H2O"]
_WILKE_ORDER = WILKE_ORDER
_WILKE_IDX = {s: i + 1 for i, s in enumerate(_WILKE_ORDER)}

# Fuller diffusion volumes [cm³/mol] (gas-phase binary estimate).
FULLER_VOL_CM3_MOL: dict[str, float] = {
    "N2": 18.5,
    "CO2": 26.9,
    "H2": 6.12,
    "CH4": 25.14,
    "H2O": 13.1,
}


def _diff_ij_key(sp_a: str, sp_b: str) -> str:
    return f"Diff{_WILKE_IDX[sp_a]}{_WILKE_IDX[sp_b]}"


def diff_ij_source(cfg: DfmConfig) -> str:
    """
    Normalize Diff_ij source: ``fuller`` | ``user``.

    Legacy aliases for user-supplied table: calibrated, table, manual, prescribed.
    """
    tr = cfg.gas_phase.get("transport") or {}
    raw = str(tr.get("diff_ij_source") or "fuller").strip().lower()
    if raw in ("user", "calibrated", "table", "manual", "prescribed", "maxwell_stefan", "ms"):
        return "user"
    if raw == "fuller":
        return "fuller"
    return raw


def fuller_binary_m2_s(
    T_k: float,
    P_pa: float,
    sp_a: str,
    sp_b: str,
    mm_g_mol: dict[str, float],
) -> float:
    """Fuller et al. gas binary diffusivity [m²/s]."""
    if sp_a not in FULLER_VOL_CM3_MOL or sp_b not in FULLER_VOL_CM3_MOL:
        raise ValueError(f"Fuller volumes missing for pair {sp_a}-{sp_b}.")
    p_atm = float(P_pa) / 101_325.0
    if p_atm <= 0:
        raise ValueError("Pressure must be positive for Fuller Diff_ij.")
    T_k = float(T_k)
    if T_k <= 0:
        raise ValueError("Temperature must be positive for Fuller Diff_ij.")
    va = FULLER_VOL_CM3_MOL[sp_a]
    vb = FULLER_VOL_CM3_MOL[sp_b]
    ma = float(mm_g_mol[sp_a])
    mb = float(mm_g_mol[sp_b])
    d_cm2_s = (
        0.00143
        * (T_k**1.75)
        / (
            p_atm
            * (va ** (1.0 / 3.0) + vb ** (1.0 / 3.0)) ** 2
            * math.sqrt(1.0 / ma + 1.0 / mb)
        )
    )
    return d_cm2_s * 1e-4


def fuller_binary_diffusion_map(
    cfg: DfmConfig,
    *,
    T_k: float | None = None,
    P_pa: float | None = None,
) -> dict[str, float]:
    """All off-diagonal Diff_ij [m²/s] from Fuller at (T, P)."""
    T = float(T_k) if T_k is not None else T_isothermal_K(cfg)
    if P_pa is not None:
        P = float(P_pa)
    else:
        P = float(cfg.gas_phase.get("P_total_Pa") or 0.0)
    if P <= 0:
        raise ValueError("gas_phase.P_total_Pa required for Fuller Diff_ij.")
    mm = {s: molar_mass_kg(s, cfg) * 1000.0 for s in _WILKE_ORDER}
    out: dict[str, float] = {}
    for i, a in enumerate(_WILKE_ORDER):
        for j, b in enumerate(_WILKE_ORDER):
            if i == j:
                continue
            out[_diff_ij_key(a, b)] = fuller_binary_m2_s(T, P, a, b, mm)
    return out


def binary_diffusion_map(cfg: DfmConfig) -> dict[str, float]:
    """
    Diff_ij [m²/s] for Wilke.

    - fuller: computed at runtime from Fuller(T_isothermal, P_total)
    - user: transport.binary_diffusion_m2_s (required)
    """
    tr = cfg.gas_phase.get("transport") or {}
    source = diff_ij_source(cfg)
    if source == "fuller":
        return fuller_binary_diffusion_map(cfg)
    if source == "user":
        per = tr.get("binary_diffusion_m2_s") or {}
        out = {str(k): float(v) for k, v in per.items() if float(v) > 0}
        if out:
            return out
        raise ValueError(
            "diff_ij_source=user requires transport.binary_diffusion_m2_s "
            "with positive Diff_ij entries (values you supply in the table)."
        )
    # Unknown: prefer table, else Fuller
    per = tr.get("binary_diffusion_m2_s") or {}
    out = {str(k): float(v) for k, v in per.items() if float(v) > 0}
    if out:
        return out
    return fuller_binary_diffusion_map(cfg)


def _lookup_diff_ij(sp_a: str, sp_b: str, diff_map: dict[str, float]) -> float:
    key = _diff_ij_key(sp_a, sp_b)
    if key in diff_map:
        return float(diff_map[key])
    rev = _diff_ij_key(sp_b, sp_a)
    if rev in diff_map:
        return float(diff_map[rev])
    raise KeyError(f"Missing binary diffusion entry for pair {sp_a}-{sp_b} ({key}).")


def c_min_mol_m3(cfg: DfmConfig) -> float:
    tr = cfg.gas_phase.get("transport") or {}
    return max(float(tr.get("c_min_mol_m3", 1e-7)), 0.0)


def smooth_concentration(c: float, c_min: float) -> float:
    """COMSOL-style regularization: 0.5*(c + sqrt(c² + cmin²))."""
    c = float(c)
    cm = float(c_min)
    if cm <= 0:
        return max(c, 0.0)
    return 0.5 * (c + (c * c + cm * cm) ** 0.5)


def mole_fractions_from_C(C: dict[str, float], cfg: DfmConfig) -> dict[str, float]:
    """Mole fractions y_i from molar concentrations (Wilke species set)."""
    c_min = c_min_mol_m3(cfg)
    active = [s for s in _WILKE_ORDER if s in C]
    if not active:
        active = list(_WILKE_ORDER)
    c_s = {s: smooth_concentration(float(C.get(s, 0.0)), c_min) for s in active}
    c_tot = sum(c_s.values())
    if c_tot <= 0:
        inert = str(cfg.gas_phase.get("inert_species", "N2"))
        return {s: 1.0 if s == inert else 0.0 for s in _WILKE_ORDER}
    return {s: c_s.get(s, 0.0) / c_tot for s in _WILKE_ORDER}


def molar_mass_g_mol(species: str, cfg: DfmConfig) -> float:
    return molar_mass_kg(species, cfg) * 1000.0


def wilke_Dm_m2_s(
    species: str,
    C: dict[str, float],
    cfg: DfmConfig,
    *,
    diff_map: dict[str, float] | None = None,
) -> float:
    """
    Wilke Dm_i [m²/s] for species in _WILKE_ORDER (matches COMSOL DmCO2, …).
    Falls back to N2 binary if species not in Wilke set.
    """
    if species not in _WILKE_ORDER:
        species = "N2"
    diff_map = diff_map or binary_diffusion_map(cfg)
    y = mole_fractions_from_C(C, cfg)
    mm = {s: molar_mass_g_mol(s, cfg) for s in _WILKE_ORDER}
    mn = sum(y[s] * mm[s] for s in _WILKE_ORDER)
    w = {s: y[s] * mm[s] / max(mn, 1e-30) for s in _WILKE_ORDER}

    denom = 0.0
    for other in _WILKE_ORDER:
        if other == species:
            continue
        d_ij = _lookup_diff_ij(species, other, diff_map)
        denom += y[other] / max(d_ij, 1e-30)
    if denom <= 0:
        per = (cfg.gas_phase.get("transport") or {}).get("D_m_per_species") or {}
        if species in per:
            return float(per[species])
        return _lookup_diff_ij(species, "CO2", diff_map)
    return (1.0 - w[species]) / denom


def make_wilke_dm_evaluator(cfg: DfmConfig, species: list[str]):
    """
    Vectorized Wilke Dm over nodal arrays (same math as wilke_Dm_m2_s).

    Returns a callable {sp: C array [mol/m³]} → {sp: Dm array [m²/s]} for the
    requested species (each must be in the Wilke set).
    """
    diff_map = binary_diffusion_map(cfg)
    c_min = c_min_mol_m3(cfg)
    mm = {s: molar_mass_g_mol(s, cfg) for s in _WILKE_ORDER}
    dij = {
        s: {o: _lookup_diff_ij(s, o, diff_map) for o in _WILKE_ORDER if o != s}
        for s in species
    }

    def evaluate(C: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        cs = {}
        for s in _WILKE_ORDER:
            c = np.asarray(C.get(s, 0.0), dtype=float)
            if c_min > 0:
                cs[s] = 0.5 * (c + np.sqrt(c * c + c_min * c_min))
            else:
                cs[s] = np.maximum(c, 0.0)
        c_tot = sum(cs.values())
        y = {s: cs[s] / c_tot for s in _WILKE_ORDER}
        mn = sum(y[s] * mm[s] for s in _WILKE_ORDER)
        out = {}
        for s in species:
            w_s = y[s] * mm[s] / np.maximum(mn, 1e-30)
            denom = sum(y[o] / max(d, 1e-30) for o, d in dij[s].items())
            out[s] = (1.0 - w_s) / denom
        return out

    return evaluate


def D_m_mode(cfg: DfmConfig) -> str:
    tr = cfg.gas_phase.get("transport") or {}
    mode = tr.get("D_m_mode", "constant")
    return str(mode).strip().lower()


def dl_dispersion_m2_s(
    D_m: float,
    dp_m: float,
    vel_d_m_s: float,
    *,
    alpha: float = 0.7,
    beta: float = 0.5,
) -> float:
    return alpha * float(D_m) + beta * float(dp_m) * max(float(vel_d_m_s), 0.0)
