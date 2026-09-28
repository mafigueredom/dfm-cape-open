"""
CAPE-OPEN DFM unit-parameter patch (step 3).

Writes geometry, CSC+LDF, methanation, cycle windows, and outlet SOPDT
onto a JSON copy. Mesh node counts stay frozen. Inlet SOPDT stays off.
"""

from __future__ import annotations

import math
from copy import deepcopy
from typing import Any

from adr.cape.contract import CapeContractError, REQUIRED_FEEDS, SOPDT_SPECIES, SPECIES

_PHASE_ROLE_ALIASES: dict[str, tuple[str, ...]] = {
    "purge2": ("purge2", "purg2"),
    "purge": ("purge", "purg"),
    "rxn": ("rxn", "reaction", "methanation"),
    "ads": ("ads", "adsorption"),
}
_CYCLE_ROLE_ORDER = ("ads", "purge", "rxn", "purge2")
_DURATION_KEYS = {
    "ads": "t_ads",
    "purge": "t_purge",
    "rxn": "t_rxn",
    "purge2": "t_purge2",
}
_DURATION_JSON = {
    "ads": "d_t_ads_s",
    "purge": "d_t_purg_s",
    "rxn": "d_t_rxn_s",
    "purge2": "d_t_purg2_s",
}
_PHASE_DEFAULTS = {
    "ads": {
        "id": "adsorption_s1",
        "adsorption_active": True,
        "methanation_active": False,
    },
    "purge": {
        "id": "purge",
        "adsorption_active": False,
        "methanation_active": False,
    },
    "rxn": {
        "id": "reaction_s2",
        "adsorption_active": False,
        "methanation_active": True,
    },
    "purge2": {
        "id": "purge2",
        "adsorption_active": False,
        "methanation_active": False,
    },
}


def classify_phase(phase: dict[str, Any]) -> str | None:
    """Map a cycle phase onto a convention-B feed id."""
    pid = str(phase.get("id") or "").lower()
    for role, aliases in _PHASE_ROLE_ALIASES.items():
        if any(a in pid for a in aliases):
            return role
    if phase.get("methanation_active"):
        return "rxn"
    if phase.get("adsorption_active"):
        return "ads"
    return None


def cycle_phase_windows(raw: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Positive-duration phases keyed by feed role (first match wins)."""
    phases = ((raw.get("simulation") or {}).get("cycle") or {}).get("phases") or []
    out: dict[str, dict[str, Any]] = {}
    for ph in phases:
        t0 = float(ph.get("t_start_s") or 0.0)
        t1 = float(ph.get("t_end_s") or 0.0)
        dt = t1 - t0
        if dt <= 1e-12:
            continue
        role = classify_phase(ph)
        if role is None or role in out:
            continue
        out[role] = {"phase": ph, "t_s": dt, "t_start_s": t0, "t_end_s": t1}
    return out


def _inert_inlet() -> dict[str, float]:
    return {sp: (1.0 if sp == "N2" else 0.0) for sp in SPECIES}


def _clone_phase(src: dict[str, Any] | None, role: str) -> dict[str, Any]:
    flags = _PHASE_DEFAULTS[role]
    base = deepcopy(src) if src else {}
    base["id"] = flags["id"]
    base["adsorption_active"] = flags["adsorption_active"]
    base["methanation_active"] = flags["methanation_active"]
    if "inlet_mole_fractions" not in base:
        base["inlet_mole_fractions"] = _inert_inlet()
    if "inlet_C_mol_m3" not in base:
        base["inlet_C_mol_m3"] = {sp: 0.0 for sp in SPECIES}
    return base


def _existing_durations(raw: dict[str, Any]) -> dict[str, float]:
    sim = raw.get("simulation") or {}
    stored = sim.get("cycle_durations") or sim.get("comsol_cycle") or {}
    windows = cycle_phase_windows(raw)
    out: dict[str, float] = {}
    for role, json_key in _DURATION_JSON.items():
        if role in windows:
            out[role] = float(windows[role]["t_s"])
        else:
            try:
                out[role] = float(stored.get(json_key) or 0.0)
            except (TypeError, ValueError):
                out[role] = 0.0
    return out


def rebuild_cycle_phases(raw: dict[str, Any], cycle: dict[str, Any]) -> dict[str, Any]:
    """Rebuild consecutive phase windows from durations. Skip a role if t=0."""
    durs = _existing_durations(raw)
    for role, key in _DURATION_KEYS.items():
        if key in cycle:
            durs[role] = float(cycle[key])
    if all(v <= 0.0 for v in durs.values()):
        raise CapeContractError("cycle durations are all zero")

    sim = raw.setdefault("simulation", {})
    cycle_block = sim.setdefault("cycle", {})
    existing = list(cycle_block.get("phases") or [])
    by_role: dict[str, dict[str, Any]] = {}
    for ph in existing:
        role = classify_phase(ph)
        if role and role not in by_role:
            by_role[role] = ph

    ordered: list[dict[str, Any]] = []
    t = 0.0
    for role in _CYCLE_ROLE_ORDER:
        dt = float(durs[role])
        if dt <= 0.0:
            continue
        src = by_role.get(role)
        if src is None and role == "purge2":
            src = by_role.get("purge")
        ph = _clone_phase(src, role)
        ph["t_start_s"] = t
        ph["t_end_s"] = t + dt
        ordered.append(ph)
        t += dt

    missing = [fid for fid in REQUIRED_FEEDS if all(classify_phase(p) != fid for p in ordered)]
    if missing:
        raise CapeContractError(
            f"cycle rebuild dropped required phase(s) {missing}; "
            "t_ads, t_purge, and t_rxn must be > 0"
        )

    n = int(cycle["n_cycles"]) if "n_cycles" in cycle else max(
        1, int(sim.get("n_cycles") or 1)
    )
    cycle_block["phases"] = ordered
    sim["n_cycles"] = n
    sim["t_end_s"] = float(n) * t
    sim["cycle_durations"] = {
        **(sim.get("cycle_durations") or {}),
        "d_t_ads_s": durs["ads"],
        "d_t_purg_s": durs["purge"],
        "d_t_rxn_s": durs["rxn"],
        "d_t_purg2_s": durs["purge2"],
    }
    sim.pop("comsol_cycle", None)
    if "adsorption_end_s" in cycle_block:
        cycle_block["adsorption_end_s"] = durs["ads"]
    steps = sim.get("comsol_steps")
    if isinstance(steps, dict) and "step1_CO2_OFF_location_s" in steps:
        steps["step1_CO2_OFF_location_s"] = durs["ads"]
    return {
        "n_cycles": n,
        "t_cycle_s": t,
        "t_horizon_s": float(n) * t,
        "durations_s": dict(durs),
        "phase_ids": [p.get("id") for p in ordered],
    }


def _apply_geometry(raw: dict[str, Any], geom: dict[str, float]) -> dict[str, Any]:
    g = raw.setdefault("geometry", {})
    bed = raw.setdefault("bed_structure", {})
    written: dict[str, float] = {}
    if "L" in geom:
        g["L_m"] = geom["L"]
        written["L"] = geom["L"]
    if "R" in geom:
        g["R_m"] = geom["R"]
        g["D_m"] = 2.0 * geom["R"]
        written["R"] = geom["R"]
        written["D"] = 2.0 * geom["R"]
    if "d_p" in geom:
        bed["d_p_m"] = geom["d_p"]
        written["d_p"] = geom["d_p"]
    if "epsilon" in geom:
        bed["epsilon"] = geom["epsilon"]
        written["epsilon"] = geom["epsilon"]
    if "rho_a" in geom:
        bed["rho_a_kg_m3_bed"] = geom["rho_a"]
        written["rho_a"] = geom["rho_a"]

    L = float(g.get("L_m") or 0.0)
    R = float(g.get("R_m") or 0.0)
    if "R" in geom:
        g["D_m"] = 2.0 * R
    mesh = raw.get("mesh") or {}
    nodes = mesh.get("nodes") or {}
    eff = mesh.setdefault("spacing_effective_m", {}) if "mesh" in raw else {}
    if "L" in geom and int(nodes.get("Nz") or 0) >= 2:
        eff["dz"] = L / (int(nodes["Nz"]) - 1)
    if "R" in geom and int(nodes.get("Nr") or 0) >= 2:
        eff["dr"] = R / (int(nodes["Nr"]) - 1)

    eps = float(bed.get("epsilon") or 0.0)
    rho_a = float(bed.get("rho_a_kg_m3_bed") or 0.0)
    if 0.0 < eps < 1.0 and rho_a > 0.0:
        bed["rho_s_kg_m3"] = rho_a / (1.0 - eps)
    if L > 0.0 and R > 0.0 and rho_a > 0.0:
        bed["m_dfm_kg"] = rho_a * math.pi * R * R * L

    notes: list[str] = []
    d_p = float(bed.get("d_p_m") or 0.0)
    D = float(g.get("D_m") or (2.0 * R))
    if d_p > 0.0 and D > 0.0 and eps > 0.0:
        corr = 0.375 + 0.34 * (d_p / D)
        if abs(eps - corr) > 1e-3:
            notes.append(
                f"epsilon={eps:.6g} is typed (not recomputed); "
                f"packing correlation 0.375+0.34(d_p/D)={corr:.6g}"
            )
    return {"written": written, "notes": notes}


def _apply_adsorption(raw: dict[str, Any], ads: dict[str, float]) -> dict[str, float]:
    block = (raw.setdefault("solid_kinetics", {})).setdefault("adsorption_CO2", {})
    csc = block.setdefault("csc", {})
    lang = block.setdefault("langmuir", {})
    written: dict[str, float] = {}
    if "k_ads" in ads:
        block["k_ads"] = ads["k_ads"]
        lang["k_ads"] = ads["k_ads"]
        written["k_ads"] = ads["k_ads"]
    if "k_des" in ads:
        block["k_des"] = ads["k_des"]
        written["k_des"] = ads["k_des"]
    if "q_e" in ads:
        csc["q_e_mol_kg"] = ads["q_e"]
        lang["q_sat_mol_kg"] = ads["q_e"]
        written["q_e"] = ads["q_e"]
    if "K_C" in ads:
        csc["K_C_Pa_inv"] = ads["K_C"]
        csc["K_c_Pa_inv"] = ads["K_C"]
        lang["K_Pa_inv"] = ads["K_C"]
        written["K_C"] = ads["K_C"]
    if "K_R" in ads:
        csc["K_R"] = ads["K_R"]
        written["K_R"] = ads["K_R"]
    if "a_csc" in ads:
        csc["a"] = ads["a_csc"]
        written["a_csc"] = ads["a_csc"]
    return written


def _apply_reaction(raw: dict[str, Any], rxn: dict[str, float]) -> dict[str, float]:
    meth = (raw.setdefault("solid_kinetics", {})).setdefault("methanation", {})
    arr = meth.setdefault("arrhenius", {})
    pl = meth.setdefault("power_law", {})
    written: dict[str, float] = {}
    if "A_rxn" in rxn:
        arr["A_rxn"] = rxn["A_rxn"]
        pl["k_ref"] = rxn["A_rxn"]
        written["A_rxn"] = rxn["A_rxn"]
    if "Ea" in rxn:
        arr["Ea_J_mol"] = rxn["Ea"]
        written["Ea"] = rxn["Ea"]
    if "m1" in rxn:
        arr["m1"] = rxn["m1"]
        pl["n_H2"] = rxn["m1"]
        written["m1"] = rxn["m1"]
    if "m2" in rxn:
        arr["m2"] = rxn["m2"]
        pl["n_CO2"] = rxn["m2"]
        written["m2"] = rxn["m2"]
    return written


def _apply_outlet_sopdt(raw: dict[str, Any], sopdt: dict[str, Any]) -> dict[str, Any]:
    z0 = (
        (raw.setdefault("gas_phase", {}))
        .setdefault("boundaries_mass", {})
        .setdefault("z_equals_0", {})
    )
    inlet = z0.setdefault("sopdt", {})
    inlet["enabled"] = False
    filt = z0.setdefault("outlet_sopdt_filter", {})
    written: dict[str, Any] = {}
    if "enabled" in sopdt:
        filt["enabled"] = bool(sopdt["enabled"])
        written["enabled"] = bool(sopdt["enabled"])
    per = sopdt.get("per_species") or {}
    dest = filt.setdefault("per_species", {})
    for sp, row in per.items():
        if sp not in SOPDT_SPECIES:
            continue
        slot = dest.setdefault(sp, {})
        written[sp] = {}
        for side in ("rising", "falling"):
            quad = row.get(side)
            if not quad:
                continue
            slot[side] = dict(quad)
            written[sp][side] = dict(quad)
        if sp == "CO2" and "rising" in row:
            for k, val in row["rising"].items():
                filt[k] = val
    return written


def apply_unit_parameters(
    raw: dict[str, Any], parameters: dict[str, Any] | None
) -> dict[str, Any]:
    """
    Patch unit parameters onto ``raw`` (in place).

    ``epsilon`` is typed. ``D_m = 2R``. Inlet ``sopdt.enabled`` is forced false.
    """
    params = parameters or {}
    info: dict[str, Any] = {"notes": []}
    z0 = (
        (raw.setdefault("gas_phase", {}))
        .setdefault("boundaries_mass", {})
        .setdefault("z_equals_0", {})
    )
    z0.setdefault("sopdt", {})["enabled"] = False

    if params.get("geometry"):
        geo = _apply_geometry(raw, params["geometry"])
        info["geometry"] = geo["written"]
        info["notes"].extend(geo["notes"])
    if params.get("adsorption"):
        info["adsorption"] = _apply_adsorption(raw, params["adsorption"])
    if params.get("reaction"):
        info["reaction"] = _apply_reaction(raw, params["reaction"])
    if "cycle" in params:
        info["cycle"] = rebuild_cycle_phases(raw, params["cycle"])
    if params.get("outlet_sopdt"):
        info["outlet_sopdt"] = _apply_outlet_sopdt(raw, params["outlet_sopdt"])
    return info
