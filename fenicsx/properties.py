"""Gas/solid properties and transport coefficients from dfm_config (isothermal v1)."""

from __future__ import annotations

import math
from typing import Any

from dfm_config import DfmConfig


def T_isothermal_K(cfg: DfmConfig) -> float:
    T = cfg.energy_thermal.get("T_isothermal_K")
    if T is None or float(T) <= 0:
        raise ValueError("energy_thermal.T_isothermal_K must be set (block D).")
    return float(T)


def species_list(cfg: DfmConfig) -> list[str]:
    sp = cfg.gas_phase.get("species")
    if not sp:
        raise ValueError("gas_phase.species missing in config.")
    return list(sp)


def inlet_bc_mode(cfg: DfmConfig) -> str:
    """
    Inlet mass BC for diluted species at z = 0.

    Returns 'danckwerts' (COMSOL FluxDanckwerts) or 'dirichlet' (default).

    Preference order: fenics_bc → inlet_bc → type (type may be
    'inlet_with_sopdt', which selects the schedule, not the BC family).
    """
    bm = cfg.gas_phase.get("boundaries_mass") or {}
    z0 = bm.get("z_equals_0") or {}
    for key in ("fenics_bc", "inlet_bc", "type"):
        raw = z0.get(key)
        if raw is None:
            continue
        mode = str(raw).strip().lower().replace("-", "_")
        if mode in (
            "flux_danckwerts",
            "danckwerts",
            "inflow_danckwerts",
            "comsol_fluxdanckwerts",
        ):
            return "danckwerts"
        if mode in ("dirichlet", "concentration"):
            return "dirichlet"
    return "dirichlet"


def R_gas(cfg: DfmConfig) -> float:
    return float(cfg.gas_phase.get("R_J_mol_K", 8.314462618))


# Molar masses [kg/mol] — override via gas_phase.molar_mass_kg_mol in config.
_DEFAULT_MOLAR_MASS_KG_MOL: dict[str, float] = {
    "CO2": 44.009e-3,
    "H2": 2.016e-3,
    "CH4": 16.043e-3,
    "H2O": 18.015e-3,
    "N2": 28.014e-3,
}


def molar_mass_kg(species: str, cfg: DfmConfig) -> float:
    per = (cfg.gas_phase.get("molar_mass_kg_mol") or {}).get(species)
    if per is not None:
        return float(per)
    if species not in _DEFAULT_MOLAR_MASS_KG_MOL:
        raise ValueError(f"No molar mass for species {species!r}.")
    return _DEFAULT_MOLAR_MASS_KG_MOL[species]


def rho_a_bed_kg_m3(cfg: DfmConfig) -> float:
    """Solid apparent density ρ_a [kg_s / m³_bed] (Darcy mass source prefactor)."""
    bed = cfg.raw.get("bed_structure") or {}
    rho_a = bed.get("rho_a_kg_m3_bed")
    if rho_a is not None and float(rho_a) > 0:
        return float(rho_a)
    eps = float(cfg.geometry.epsilon)
    rho = float(cfg.geometry.rho_s_kg_m3)
    return (1.0 - eps) * rho


def rho_a_over_epsilon(cfg: DfmConfig) -> float:
    """ρ_a/ε = (1−ε) ρ_s / ε [kg_s / m³_void] — diluted-species molar source prefactor."""
    eps = float(cfg.geometry.epsilon)
    return rho_a_bed_kg_m3(cfg) / max(eps, 1e-12)


def rho_a_solid_bulk(cfg: DfmConfig) -> float:
    """Alias for ρ_a/ε (diluted-species transport)."""
    return rho_a_over_epsilon(cfg)


def fluid_density(C: dict[str, float], cfg: DfmConfig) -> float:
    """Mixture mass density ρ_f = Σ_i C_i M_i [kg/m³] from molar concentrations."""
    return sum(float(C[sp]) * molar_mass_kg(sp, cfg) for sp in species_list(cfg))


def darcy_permeability_m2(cfg: DfmConfig) -> float:
    """Permeability 𝒦 [m²] (`gas_phase.momentum.permeability_m2`)."""
    mom = cfg.gas_phase.get("momentum") or {}
    K = mom.get("permeability_m2")
    if K is None or float(K) <= 0:
        raise ValueError("gas_phase.momentum.permeability_m2 must be positive for Darcy law.")
    return float(K)


def fluid_viscosity_Pa_s(cfg: DfmConfig) -> float:
    """Dynamic viscosity μ_f [Pa·s] (`gas_phase.momentum.viscosity_Pa_s` or mixture estimate)."""
    mom = cfg.gas_phase.get("momentum") or {}
    mu = mom.get("viscosity_Pa_s")
    if mu is not None and float(mu) > 0:
        return float(mu)
    # Fallback: H2-like order of magnitude at ~573 K (override in config for accuracy).
    return 1.2e-5


def total_pressure_Pa(cfg: DfmConfig) -> float:
    P = cfg.gas_phase.get("P_total_Pa")
    if P is None or float(P) <= 0:
        raise ValueError("gas_phase.P_total_Pa required for Darcy pressure profile.")
    return float(P)


def ads_feed_P_CO2_Pa(cfg: DfmConfig) -> float:
    """
    Reference CO₂ partial pressure P_CO2,0 for θ normalization [Pa].

    Prefer ads-phase inlet mole fraction: P_CO2,0 = P_total · y_CO2.
    Else use ads-phase / t=0 inlet concentration: P = C_CO2 · R · T.
    """
    P = total_pressure_Pa(cfg)
    sim = cfg.raw.get("simulation") or {}
    phases = (sim.get("cycle") or {}).get("phases") or []
    ads_phases = [ph for ph in phases if ph.get("adsorption_active")]
    candidates = ads_phases or list(phases)
    for ph in candidates:
        y = ph.get("inlet_mole_fractions") or {}
        if "CO2" in y and float(y["CO2"]) > 0:
            return P * float(y["CO2"])
        cin = ph.get("inlet_C_mol_m3") or {}
        if "CO2" in cin and float(cin["CO2"]) > 0:
            return partial_pressure(float(cin["CO2"]), cfg)
    try:
        cin0 = inlet_concentrations(cfg)
        c_co2 = float(cin0.get("CO2", 0.0))
        if c_co2 > 0:
            return partial_pressure(c_co2, cfg)
    except ValueError:
        pass
    raise ValueError(
        "Cannot form P_CO2,0 for q_e,const = CSC(P_CO2,0): need ads-phase "
        "inlet y_CO2 > 0 or C_CO2 > 0."
    )


def q_e_const_mol_kg(cfg: DfmConfig) -> float:
    """
    Coverage normalizer q_{e,const} for θ = q / q_{e,const}.

    Always COMSOL-style qe_max = CSC(P_CO2,0). Stored JSON
    ``qe_const_mol_kg`` is ignored (legacy; recomputed from CSC + ads feed).
    """
    # Lazy import: kinetics imports this module at load time.
    from kinetics import q_e_csc

    return max(float(q_e_csc(ads_feed_P_CO2_Pa(cfg), cfg)), 1e-15)


def q_e_max_mol_kg(cfg: DfmConfig) -> float:
    """CSC capacity scale q_e [mol/kg] (not the θ normalizer)."""
    ads = cfg.solid_kinetics.get("adsorption_CO2") or {}
    csc = ads.get("csc") or {}
    q_e = csc.get("q_e_mol_kg")
    if q_e is None or float(q_e) <= 0:
        raise ValueError("adsorption_CO2.csc.q_e_mol_kg > 0 is required.")
    return float(q_e)


def superficial_velocity(cfg: DfmConfig, Q_m3_s: float | None = None) -> float:
    """Superficial Darcy velocity Vel_D = Q / (pi R^2) [m/s] (COMSOL dl.w)."""
    vel = cfg.gas_phase.get("velocity") or {}
    if Q_m3_s is not None and float(Q_m3_s) > 0:
        Q = float(Q_m3_s)
    else:
        Q = vel.get("Q_m3_s_at_T_ref")
        if Q is None or float(Q) <= 0:
            raise ValueError("velocity.Q_m3_s_at_T_ref required for superficial velocity.")
        Q = float(Q)
    R = cfg.geometry.R_m
    return Q / (math.pi * R * R)


def D_m_species(species: str, cfg: DfmConfig, C: dict[str, float] | None = None) -> float:
    """Molecular diffusivity D_m [m²/s] (constant config or Wilke from local C)."""
    from fenicsx.wilke_transport import D_m_mode, wilke_Dm_m2_s

    if D_m_mode(cfg) == "wilke" and C is not None:
        return wilke_Dm_m2_s(species, C, cfg)
    tr = cfg.gas_phase.get("transport") or {}
    per = tr.get("D_m_per_species") or {}
    if species in per and float(per[species]) > 0:
        return float(per[species])
    diff = cfg.gas_phase.get("diffusion") or {}
    row = (diff.get("D_mol_per_species") or {}).get(species)
    if row:
        D_ref = row.get("D_mol_ref_m2_s")
        T_ref = row.get("T_ref_K")
        n = row.get("T_exponent_n")
        if D_ref is not None and T_ref is not None and n is not None:
            T = T_isothermal_K(cfg)
            return float(D_ref) * (T / float(T_ref)) ** float(n)
    # Fallback: undo tortuosity/porosity scaling from effective D if only D_eff is calibrated.
    D_eff = D_eff_species(species, cfg)
    eps = float(cfg.geometry.epsilon)
    tau = float(diff.get("tau") or 1.0)
    return D_eff * float(tau) / max(eps, 1e-12)


def D_z_species(
    species: str,
    u_magnitude: float,
    cfg: DfmConfig,
    *,
    C: dict[str, float] | None = None,
    vel_d_m_s: float | None = None,
) -> float:
    """
    Axial dispersion D_z [m²/s].

    Modes (gas_phase.transport.D_z_mode):
      dispersion_formula — alpha D_m + beta d_p |u| (COMSOL DL(); default alpha=0.7, beta=0.5)
      constant           — fixed D_z per species (legacy shortcut; not POW_V07 active model)
    """
    tr = cfg.gas_phase.get("transport") or {}
    mode = str(tr.get("D_z_mode", "dispersion_formula")).strip()
    if mode == "constant":
        per = tr.get("D_z_const_per_species") or {}
        if species in per and float(per[species]) > 0:
            return float(per[species])
        fallback = tr.get("D_z_const_m2_s")
        if fallback is not None and float(fallback) > 0:
            return float(fallback)
        raise ValueError(f"transport.D_z_const_per_species[{species!r}] required for constant mode.")
    return D_z_dispersion(species, u_magnitude, cfg, C=C, vel_d_m_s=vel_d_m_s)


def D_z_dispersion(
    species: str,
    u_magnitude: float,
    cfg: DfmConfig,
    *,
    C: dict[str, float] | None = None,
    vel_d_m_s: float | None = None,
) -> float:
    """
    D_z = alpha D_m + beta d_p |u| [m²/s].

    COMSOL POW_V07: D_c_i = DL(Dm_i, dp_dfm, dl.w) with dl.w = superficial Darcy velocity.
    Pass vel_d_m_s directly when known; otherwise derive from interstitial |u| and basis flag.
    """
    from fenicsx.wilke_transport import dl_dispersion_m2_s

    tr = cfg.gas_phase.get("transport") or {}
    alpha = float(tr.get("dispersion_alpha", 0.7))
    beta = float(tr.get("dispersion_beta", 0.5))
    d_p = cfg.geometry.d_p_m
    if d_p is None or float(d_p) <= 0:
        raise ValueError("bed_structure.d_p_m required for D_z = alpha*D_m + beta*d_p*|u|.")
    if vel_d_m_s is not None:
        u_superficial = max(float(vel_d_m_s), 0.0)
    else:
        basis = str(tr.get("D_z_velocity_basis", "interstitial")).strip().lower()
        u_superficial = float(u_magnitude)
        if basis not in ("superficial", "darcy", "dl.w", "vel_d"):
            u_superficial = u_superficial * float(cfg.geometry.epsilon)
    D_m = D_m_species(species, cfg, C=C)
    return dl_dispersion_m2_s(D_m, float(d_p), u_superficial, alpha=alpha, beta=beta)


def diluted_species_initial_C(cfg: DfmConfig, species: list[str]) -> dict[str, float]:
    """
    Initial concentrations from model sheet:
        P = P_op,  C_N2 = P_op / (R T_op),  other species 0,  θ_CO2 = 0.
    """
    sim = cfg.raw.get("simulation") or {}
    ic = sim.get("initial_C_mol_m3")
    if ic:
        return {sp: float(ic.get(sp, 0.0)) for sp in species}
    P = total_pressure_Pa(cfg)
    R = R_gas(cfg)
    T = T_isothermal_K(cfg)
    c_inert = P / (R * T)
    inert = str(cfg.gas_phase.get("inert_species", "N2"))
    out = {sp: 0.0 for sp in species}
    if inert in out:
        out[inert] = c_inert
    return out


def partial_pressure(C_mol_m3: float, cfg: DfmConfig) -> float:
    return C_mol_m3 * R_gas(cfg) * T_isothermal_K(cfg)


def linear_property(a0: float | None, a1: float | None, T: float) -> float:
    if a0 is None or a1 is None:
        raise ValueError("Property linear coefficients a0, a1 required.")
    return float(a0) + float(a1) * T


def interstitial_velocity(cfg: DfmConfig, Q_m3_s: float | None = None) -> float:
    """
    Interstitial velocity v_z = Q / (ε π R²).

    Always derived from volumetric flow Q (and bed ε, R). Optional legacy
    JSON ``v_z_m_s`` is ignored if present.
    """
    vel = cfg.gas_phase.get("velocity") or {}
    if Q_m3_s is not None and float(Q_m3_s) > 0:
        Q = float(Q_m3_s)
    else:
        Q = vel.get("Q_m3_s_at_T_ref")
        if Q is None or float(Q) <= 0:
            raise ValueError(
                "gas_phase.velocity.Q_m3_s_at_T_ref required "
                "(v_z is calculated as Q/(ε π R²); do not set v_z_m_s as input)."
            )
        Q = float(Q)
    eps = float(cfg.geometry.epsilon)
    if eps <= 0:
        raise ValueError("bed_structure.epsilon must be positive to compute v_z.")
    R = float(cfg.geometry.R_m)
    A = math.pi * R * R
    return Q / (A * eps)


def D_eff_species(species: str, cfg: DfmConfig) -> float:
    diff = cfg.gas_phase.get("diffusion") or {}
    if diff.get("use_wakao_smith"):
        raise NotImplementedError(
            "Wakao–Smith D_eff not implemented in fenicsx v1; set diffusion.mode=single "
            "with D_eff_single_m2_s or per-species D_mol and tau."
        )
    mode = diff.get("mode", "single")
    eps = cfg.geometry.epsilon
    if mode == "single":
        D = diff.get("D_eff_single_m2_s")
        if D is None or float(D) <= 0:
            raise ValueError("diffusion.D_eff_single_m2_s required for mode=single.")
        scale = float((diff.get("D_eff_scale_per_species") or {}).get(species, 1.0))
        return float(D) * scale
    per = diff.get("D_mol_per_species") or {}
    row = per.get(species)
    if not row:
        raise ValueError(f"Missing D_mol_per_species for {species}.")
    D_ref = row.get("D_mol_ref_m2_s")
    T_ref = row.get("T_ref_K")
    n = row.get("T_exponent_n")
    tau = diff.get("tau")
    if D_ref is None or T_ref is None or n is None or tau is None:
        raise ValueError(f"Incomplete D_mol/tau for {species}.")
    T = T_isothermal_K(cfg)
    D_mol = float(D_ref) * (T / float(T_ref)) ** float(n)
    return D_mol * eps / float(tau)


def inlet_concentrations(cfg: DfmConfig) -> dict[str, float]:
    sim = cfg.raw.get("simulation") or {}
    cycle = sim.get("cycle")
    if cycle and cycle.get("phases"):
        from fenicsx.cycle_scheduler import CycleScheduler

        return CycleScheduler(cfg).inlet_at(0.0)
    cin = sim.get("inlet_C_mol_m3")
    if not cin:
        raise ValueError(
            "Add simulation.inlet_C_mol_m3 to dfm_config.json "
            "(map species name -> mol/m^3) before running FEniCSx."
        )
    out = {str(k): float(v) for k, v in cin.items()}
    for sp in species_list(cfg):
        if sp not in out:
            raise ValueError(f"simulation.inlet_C_mol_m3 missing species {sp!r}.")
    return out


def time_stepping(cfg: DfmConfig) -> tuple[float, float]:
    from fenicsx.cycle_scheduler import resolved_t_end_s

    sim = cfg.raw.get("simulation") or {}
    dt = sim.get("dt_s")
    t_end = resolved_t_end_s(sim)
    if dt is None or float(dt) <= 0 or t_end <= 0:
        raise ValueError(
            "simulation.dt_s must be positive, and t_end_s or n_cycles × cycle period must be positive."
        )
    return float(t_end), float(dt)
