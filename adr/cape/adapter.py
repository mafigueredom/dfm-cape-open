"""
CAPE-OPEN DFM engine adapter (steps 2–3).

Patches unit parameters, maps convention-B port feeds onto a DFM JSON
copy, runs ``SimulationController(factory_id="methanation_dfm")``, and
builds ``cape_result_v1``. Does not take packed-bed ``D_z`` from a
property package and does not write a stream temperature into
``T_isothermal_K``.
"""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path
from typing import Any

from adr.cape.contract import (
    ALLOWED_FACTORY_ID,
    REQUEST_SCHEMA,
    REQUIRED_FEEDS,
    SPECIES,
    CapeContractError,
    build_cape_result,
    parse_cape_request,
    validate_factory_id,
)
from adr.cape.balance import force_dfm_balance_enabled
from adr.cape.parameters import apply_unit_parameters, cycle_phase_windows
from dfm_config import DfmConfig, dfm_config_from_dict

_V1_ROOT = Path(__file__).resolve().parents[2]
_R_DEFAULT = 8.314462618
_P_EQUAL_TOL_PA = 1.0
_Q_WARN_RATIO = 1.25


def _as_parsed_request(request: dict[str, Any]) -> dict[str, Any]:
    if request.get("schema") == REQUEST_SCHEMA and "feeds" in request:
        if "factory_id" in request:
            validate_factory_id(request["factory_id"])
        if request.get("feeds") and isinstance(
            next(iter(request["feeds"].values()), None), dict
        ):
            first = next(iter(request["feeds"].values()))
            if "y" in first and "CO2" in first["y"]:
                return request
    return parse_cape_request(request)


def _resolve_config_path(path: str, base_dir: Path) -> Path:
    p = Path(path)
    if p.is_file():
        return p
    cand = base_dir / path
    if cand.is_file():
        return cand
    raise CapeContractError(f"DFM config not found: {path}")


def _is_fermentation_config(data: dict[str, Any]) -> bool:
    """DFM-only check; does not import ``adr.fermentation``."""
    schema = str(data.get("schema") or data.get("schema_version") or "")
    if schema.lower().startswith("fermentation"):
        return True
    meta = data.get("meta") or {}
    sim = data.get("simulation") or {}
    fid = str(sim.get("factory_id") or meta.get("factory_id") or "")
    return fid == "fermentation"


def _load_base_raw(request: dict[str, Any], base_dir: Path) -> dict[str, Any]:
    inline = request.get("config")
    if isinstance(inline, dict) and "geometry" in inline and "simulation" in inline:
        data = deepcopy(inline)
    else:
        rel = request.get("config_path")
        if not rel:
            raise CapeContractError("provide config_path or inline DFM config")
        path = _resolve_config_path(str(rel), base_dir)
        data = json.loads(path.read_text(encoding="utf-8"))
    if _is_fermentation_config(data):
        raise CapeContractError(
            "config is fermentation; CAPE adapter accepts only methanation_dfm"
        )
    fid = str((data.get("simulation") or {}).get("factory_id") or "")
    if fid and fid != ALLOWED_FACTORY_ID:
        validate_factory_id(fid)
    return data


def cycle_period_s(windows: dict[str, dict[str, Any]]) -> float:
    if not windows:
        raise CapeContractError("simulation.cycle.phases has no positive-duration steps")
    t0 = min(float(w["t_start_s"]) for w in windows.values())
    t1 = max(float(w["t_end_s"]) for w in windows.values())
    period = t1 - t0
    if period <= 0.0:
        raise CapeContractError("cycle period T must be > 0")
    return period


def bed_temperature_K(raw: dict[str, Any]) -> float:
    T = (raw.get("energy_thermal") or {}).get("T_isothermal_K")
    try:
        t = float(T)
    except (TypeError, ValueError) as exc:
        raise CapeContractError(
            "energy_thermal.T_isothermal_K is required; stream T is not used"
        ) from exc
    if t <= 0.0:
        raise CapeContractError("energy_thermal.T_isothermal_K must be > 0")
    return t


def bed_pressure_Pa(request: dict[str, Any]) -> tuple[float, list[str]]:
    """One bed P from Feed_ads. Extra port pressures are not stacked."""
    feeds = request.get("feeds") or {}
    ads = feeds.get("ads")
    if not isinstance(ads, dict):
        raise CapeContractError("feeds.ads is required for bed P")
    p_bed = float(ads["P_Pa"])
    notes: list[str] = []
    for fid, feed in feeds.items():
        if fid == "ads" or not isinstance(feed, dict):
            continue
        p_k = float(feed["P_Pa"])
        if abs(p_k - p_bed) > _P_EQUAL_TOL_PA:
            notes.append(
                f"feeds.{fid}.P_Pa={p_k:g} differs from Feed_ads P={p_bed:g}; "
                "using Feed_ads (no hydraulic stack)"
            )
    return p_bed, notes


def gas_constant(raw: dict[str, Any]) -> float:
    try:
        return float((raw.get("gas_phase") or {}).get("R_J_mol_K") or _R_DEFAULT)
    except (TypeError, ValueError):
        return _R_DEFAULT


def concentrations_from_y(
    y: dict[str, float], *, P_Pa: float, T_K: float, R: float
) -> dict[str, float]:
    ct = P_Pa / (R * T_K)
    return {sp: float(y.get(sp, 0.0)) * ct for sp in SPECIES}


def convention_b_Q(
    F_mol_s: float, *, T_s: float, t_k: float, T_K: float, P_Pa: float, R: float
) -> tuple[float, float]:
    """Return ``(F_k^bed, Q_k)``. Does not form ``T/t_k`` when ``t_k=0``."""
    if t_k <= 0.0:
        raise CapeContractError("cannot form T/t_k for a zero-duration phase")
    f_bed = float(F_mol_s) * float(T_s) / float(t_k)
    q = f_bed * R * T_K / P_Pa
    if q < 0.0:
        raise CapeContractError("convention-B Q_k must be >= 0")
    return f_bed, q


def apply_feeds_to_raw(
    raw: dict[str, Any], request: dict[str, Any]
) -> dict[str, Any]:
    """
    Patch cycle inlets from convention-B feeds.

    Writes ``inlet_C_mol_m3``, ``inlet_mole_fractions``, and ``Q_m3_s`` on
    matching phases. Sets ``gas_phase.P_total_Pa`` from Feed_ads.
    Leaves ``T_isothermal_K`` and ``gas_phase.diffusion`` unchanged.
    """
    feeds = request.get("feeds") or {}
    windows = cycle_phase_windows(raw)
    for fid in REQUIRED_FEEDS:
        if fid not in feeds:
            raise CapeContractError(f"feeds.{fid} is required")
        if fid not in windows:
            raise CapeContractError(
                f"feeds.{fid} has no matching cycle phase in the DFM config"
            )
    if "purge2" in windows and "purge2" not in feeds:
        raise CapeContractError(
            "t_purge2 > 0 requires feeds.purge2 (port connected when the phase exists)"
        )
    skipped_purge2 = False
    if "purge2" in feeds and "purge2" not in windows:
        feeds = {k: v for k, v in feeds.items() if k != "purge2"}
        skipped_purge2 = True

    T_s = cycle_period_s(windows)
    T_K = bed_temperature_K(raw)
    R = gas_constant(raw)
    P_Pa, p_notes = bed_pressure_Pa(request)

    gp = raw.setdefault("gas_phase", {})
    gp["P_total_Pa"] = P_Pa

    applied: dict[str, Any] = {}
    q_ref: float | None = None
    for fid, feed in feeds.items():
        win = windows[fid]
        t_k = float(win["t_s"])
        y = {sp: float(feed["y"].get(sp, 0.0)) for sp in SPECIES}
        cin = concentrations_from_y(y, P_Pa=P_Pa, T_K=T_K, R=R)
        f_bed, q = convention_b_Q(
            float(feed["F_mol_s"]), T_s=T_s, t_k=t_k, T_K=T_K, P_Pa=P_Pa, R=R
        )
        ph = win["phase"]
        ph["inlet_mole_fractions"] = y
        ph["inlet_C_mol_m3"] = cin
        ph["Q_m3_s"] = q
        if q_ref is None:
            q_ref = q
        applied[fid] = {
            "phase_id": ph.get("id"),
            "t_s": t_k,
            "F_mol_s": float(feed["F_mol_s"]),
            "F_bed_mol_s": f_bed,
            "Q_m3_s": q,
            "y": y,
            "inlet_C_mol_m3": cin,
        }

    vel = gp.setdefault("velocity", {})
    if q_ref is not None:
        vel["Q_m3_s_at_T_ref"] = q_ref
        vel["Fvol_L_min"] = q_ref * 1000.0 * 60.0

    notes = list(p_notes)
    if skipped_purge2:
        notes.append("feeds.purge2 ignored (t_purge2=0; port disconnected)")
    qs = [float(row["Q_m3_s"]) for row in applied.values()]
    if qs:
        q_min, q_max = min(qs), max(qs)
        if q_min > 0.0 and q_max / q_min > _Q_WARN_RATIO:
            notes.append(
                "per-phase Q differs by more than 25%; "
                "if every port carries full F_bed, velocities are T/t_k times the lab"
            )

    return {
        "T_bed_K": T_K,
        "P_Pa": P_Pa,
        "R_J_mol_K": R,
        "t_cycle_s": T_s,
        "phases": applied,
        "notes": notes,
        "note": (
            "T_isothermal_K and packed-bed D_z were not taken from the PME. "
            "Q_k = F_k*(T/t_k)*R*T_bed/P; C_i,k = y_i,k*P/(R*T_bed)."
        ),
    }


def patch_cape_request(
    request: dict[str, Any],
    *,
    base_dir: Path | str | None = None,
) -> tuple[DfmConfig, dict[str, Any]]:
    """Load the frozen DFM JSON, patch feeds, return ``(cfg, adapter_info)``."""
    req = _as_parsed_request(request)
    validate_factory_id(req.get("factory_id"))
    root = Path(base_dir) if base_dir is not None else _V1_ROOT
    raw = _load_base_raw(req, root)
    param_info = apply_unit_parameters(raw, req.get("parameters"))
    force_dfm_balance_enabled(raw)
    info = apply_feeds_to_raw(raw, req)
    info["parameters"] = param_info
    if param_info.get("notes"):
        info.setdefault("notes", []).extend(param_info["notes"])
    return dfm_config_from_dict(raw), info


def from_cape_request(
    request: dict[str, Any],
    *,
    base_dir: Path | str | None = None,
) -> DfmConfig:
    """Load the frozen DFM JSON, patch feeds, return an in-memory ``DfmConfig``."""
    cfg, _info = patch_cape_request(request, base_dir=base_dir)
    return cfg


def _attach_sopdt_series(
    result: dict[str, Any],
    run: dict[str, Any],
    cfg: DfmConfig,
) -> None:
    """Report-only SOPDT ``C(t)``. Never written onto the product stream."""
    c_t = str(result.get("C_t_export") or "both")
    if c_t not in ("sopdt", "both"):
        return
    times = run.get("times_s") or []
    raw_c = run.get("outlet_C_mol_m3") or {}
    if not times or not raw_c:
        result["outlet_C_sopdt"] = None
        result["outlet_C_sopdt_note"] = "no outlet_C_mol_m3 on the run export"
        return
    import numpy as np

    from fenicsx.inlet_sopdt import filter_outlet_series

    t = np.asarray(times, dtype=float)
    filtered: dict[str, list[float]] = {}
    for sp in SPECIES:
        if sp not in raw_c:
            continue
        y = filter_outlet_series(cfg, sp, t, np.asarray(raw_c[sp], dtype=float))
        filtered[sp] = [float(v) for v in y]
    result["outlet_C_sopdt"] = {
        "times_s": [float(x) for x in t],
        "C_mol_m3": filtered,
    }
    result["outlet_C_sopdt_note"] = (
        "SOPDT C(t) is a report series from outlet_sopdt_filter. "
        "Not used for product F_ss."
    )


def to_cape_result(
    run: dict[str, Any],
    request: dict[str, Any] | None = None,
    *,
    cfg: DfmConfig | None = None,
    adapter: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """``cape_result_v1`` from a DFM export. Product stays raw."""
    req = request or {"factory_id": ALLOWED_FACTORY_ID, "C_t_export": "both"}
    T_K = None
    P_Pa = None
    if cfg is not None:
        T_K = bed_temperature_K(cfg.raw)
        P_Pa = float((cfg.raw.get("gas_phase") or {}).get("P_total_Pa") or 101325.0)
    elif adapter is not None:
        T_K = float(adapter.get("T_bed_K") or 593.15)
        P_Pa = float(adapter.get("P_Pa") or 101325.0)
    result = build_cape_result(run, req, T_K=T_K, P_Pa=P_Pa)
    if cfg is not None:
        _attach_sopdt_series(result, run, cfg)
        from adr.cape.reports import attach_reports

        attach_reports(result)
    if adapter is not None:
        result["adapter"] = adapter
    return result


def run_cape_request(
    request: dict[str, Any],
    *,
    base_dir: Path | str | None = None,
    t_end: float | None = None,
    dt: float | None = None,
    export_json: Path | str | None = None,
    plot: bool = False,
    print_summary: bool = False,
    no_balance: bool = False,
    reports_dir: Path | str | None = None,
) -> dict[str, Any]:
    """Patch → ``SimulationController.run()`` → ``cape_result_v1``."""
    req = _as_parsed_request(request)
    cfg, adapter = patch_cape_request(req, base_dir=base_dir)

    from adr.controller import SimulationController
    from fenicsx.run_export import reactor_state_to_dict

    sim = SimulationController(factory_id=ALLOWED_FACTORY_ID).run(
        cfg,
        factory_id=ALLOWED_FACTORY_ID,
        t_end=t_end,
        dt=dt,
        export_json=export_json,
        plot=plot,
        print_summary=print_summary,
        no_balance=no_balance,
    )
    meta = {
        "fenics_implementation": (cfg.raw.get("meta") or {}).get("fenics_implementation"),
        "alignment_target": (cfg.raw.get("meta") or {}).get("alignment_target"),
        "factory_id": ALLOWED_FACTORY_ID,
    }
    run = reactor_state_to_dict(
        sim.state,
        species=sim.species,
        config_path=req.get("config_path"),
        meta=meta,
    )
    result = to_cape_result(run, req, cfg=cfg, adapter=adapter)
    if reports_dir is not None:
        from adr.cape.reports import write_cape_reports

        write_cape_reports(result, reports_dir)
    return result
