"""
CAPE-OPEN DFM JSON contract (SI).

Step 1 only: parse/validate ``cape_request_v1`` and build ``cape_result_v1``
from an existing DFM run export. Does not patch config or call FEniCS.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ALLOWED_FACTORY_ID = "methanation_dfm"
REQUEST_SCHEMA = "cape_request_v1"
RESULT_SCHEMA = "cape_result_v1"
C_T_EXPORT_MODES = ("raw", "sopdt", "both")
FEED_IDS = ("ads", "purge", "rxn", "purge2")
REQUIRED_FEEDS = ("ads", "purge", "rxn")
SPECIES = ("CO2", "H2", "CH4", "H2O", "N2")
SOPDT_SPECIES = ("CO2", "H2", "CH4", "H2O")
CAS_TO_SPECIES = {
    "124-38-9": "CO2",
    "1333-74-0": "H2",
    "74-82-8": "CH4",
    "7732-18-5": "H2O",
    "7727-37-9": "N2",
}
SPECIES_TO_CAS = {sp: cas for cas, sp in CAS_TO_SPECIES.items()}
BALANCE_TOL_DEFAULT = 1e-3
Y_SUM_TOL = 1e-4


class CapeContractError(ValueError):
    """Invalid CAPE request/result or non-DFM factory."""


def validate_factory_id(factory_id: Any) -> str:
    fid = str(factory_id or "").strip()
    if fid != ALLOWED_FACTORY_ID:
        raise CapeContractError(
            f"CAPE-OPEN DFM unit accepts only factory_id={ALLOWED_FACTORY_ID!r}, "
            f"got {fid!r}"
        )
    return fid


def _as_float(value: Any, name: str, *, positive: bool = False) -> float:
    try:
        x = float(value)
    except (TypeError, ValueError) as exc:
        raise CapeContractError(f"{name} must be a number") from exc
    if positive and x <= 0.0:
        raise CapeContractError(f"{name} must be > 0")
    if x < 0.0:
        raise CapeContractError(f"{name} must be >= 0")
    return x


def normalize_mole_fractions(y_raw: Any, *, feed_id: str) -> dict[str, float]:
    if not isinstance(y_raw, dict) or not y_raw:
        raise CapeContractError(f"feeds.{feed_id}.y must be a non-empty object")
    out = {sp: 0.0 for sp in SPECIES}
    unknown: list[str] = []
    for key, val in y_raw.items():
        k = str(key).strip()
        sp = CAS_TO_SPECIES.get(k) or (k if k in SPECIES else None)
        if sp is None:
            unknown.append(k)
            continue
        out[sp] = _as_float(val, f"feeds.{feed_id}.y[{key}]")
    if unknown:
        raise CapeContractError(
            f"feeds.{feed_id}.y has unmapped keys {unknown}; "
            f"use CAS {list(CAS_TO_SPECIES)} or {list(SPECIES)}"
        )
    total = sum(out.values())
    if abs(total - 1.0) > Y_SUM_TOL:
        raise CapeContractError(
            f"feeds.{feed_id}.y must sum to 1 (±{Y_SUM_TOL}), got {total}"
        )
    return out


def _parse_feed(raw: Any, feed_id: str, *, required: bool) -> dict[str, Any] | None:
    if raw is None:
        if required:
            raise CapeContractError(f"feeds.{feed_id} is required")
        return None
    if not isinstance(raw, dict):
        raise CapeContractError(f"feeds.{feed_id} must be an object")
    return {
        "id": feed_id,
        "F_mol_s": _as_float(raw.get("F_mol_s"), f"feeds.{feed_id}.F_mol_s"),
        "P_Pa": _as_float(raw.get("P_Pa"), f"feeds.{feed_id}.P_Pa", positive=True),
        "y": normalize_mole_fractions(raw.get("y"), feed_id=feed_id),
    }


def _parse_sopdt_quad(raw: Any, name: str) -> dict[str, float] | None:
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise CapeContractError(f"{name} must be an object")
    return {
        "Kp": _as_float(raw.get("Kp", 1.0), f"{name}.Kp", positive=True),
        "tau_s": _as_float(raw.get("tau_s", 1.0), f"{name}.tau_s", positive=True),
        "zeta": _as_float(raw.get("zeta", 1.0), f"{name}.zeta", positive=True),
        "theta_p": _as_float(raw.get("theta_p", 0.0), f"{name}.theta_p"),
    }


def _parse_optional_group(
    raw: Any,
    name: str,
    fields: dict[str, dict[str, Any]],
) -> dict[str, float]:
    if not isinstance(raw, dict):
        raise CapeContractError(f"parameters.{name} must be an object")
    out: dict[str, float] = {}
    for key, spec in fields.items():
        if key not in raw:
            continue
        out[key] = _as_float(
            raw[key], f"parameters.{name}.{key}", positive=bool(spec.get("positive"))
        )
        lo = spec.get("gt")
        hi = spec.get("lt")
        if lo is not None and out[key] <= float(lo):
            raise CapeContractError(f"parameters.{name}.{key} must be > {lo}")
        if hi is not None and out[key] >= float(hi):
            raise CapeContractError(f"parameters.{name}.{key} must be < {hi}")
    return out


def parse_unit_parameters(raw: Any) -> dict[str, Any]:
    """Validate CAPE unit parameters (geometry / kinetics / cycle / outlet SOPDT)."""
    if raw is None:
        return {}
    if not isinstance(raw, dict):
        raise CapeContractError("parameters must be an object")
    out: dict[str, Any] = {}
    if "geometry" in raw:
        out["geometry"] = _parse_optional_group(
            raw["geometry"],
            "geometry",
            {
                "L": {"positive": True},
                "R": {"positive": True},
                "d_p": {"positive": True},
                "epsilon": {"gt": 0.0, "lt": 1.0},
                "rho_a": {"positive": True},
            },
        )
    if "adsorption" in raw:
        out["adsorption"] = _parse_optional_group(
            raw["adsorption"],
            "adsorption",
            {
                "k_ads": {"positive": True},
                "k_des": {},
                "q_e": {"positive": True},
                "K_C": {"positive": True},
                "K_R": {"positive": True},
                "a_csc": {"positive": True},
            },
        )
    if "reaction" in raw:
        out["reaction"] = _parse_optional_group(
            raw["reaction"],
            "reaction",
            {
                "A_rxn": {"positive": True},
                "Ea": {},
                "m1": {"positive": True},
                "m2": {"positive": True},
            },
        )
    if "cycle" in raw:
        cyc_raw = raw["cycle"]
        if not isinstance(cyc_raw, dict):
            raise CapeContractError("parameters.cycle must be an object")
        cyc: dict[str, Any] = {}
        if "n_cycles" in cyc_raw:
            n = _as_float(cyc_raw["n_cycles"], "parameters.cycle.n_cycles", positive=True)
            n_int = int(round(n))
            if n_int < 1:
                raise CapeContractError("parameters.cycle.n_cycles must be >= 1")
            cyc["n_cycles"] = n_int
        for key in ("t_ads", "t_purge", "t_rxn", "t_purge2"):
            if key in cyc_raw:
                cyc[key] = _as_float(cyc_raw[key], f"parameters.cycle.{key}")
        out["cycle"] = cyc
    if "outlet_sopdt" in raw:
        sopdt = raw["outlet_sopdt"]
        if not isinstance(sopdt, dict):
            raise CapeContractError("parameters.outlet_sopdt must be an object")
        parsed: dict[str, Any] = {}
        if "enabled" in sopdt:
            parsed["enabled"] = bool(sopdt["enabled"])
        per = sopdt.get("per_species") or {}
        if per and not isinstance(per, dict):
            raise CapeContractError("outlet_sopdt.per_species must be an object")
        per_out: dict[str, Any] = {}
        for sp, row in (per or {}).items():
            if sp == "N2":
                raise CapeContractError("outlet_sopdt has no N2 row")
            if sp not in SOPDT_SPECIES:
                raise CapeContractError(
                    f"outlet_sopdt species {sp!r} is not a filtered DFM species"
                )
            if not isinstance(row, dict):
                raise CapeContractError(f"outlet_sopdt.per_species.{sp} must be an object")
            rising = _parse_sopdt_quad(row.get("rising"), f"outlet_sopdt.{sp}.rising")
            falling = _parse_sopdt_quad(row.get("falling"), f"outlet_sopdt.{sp}.falling")
            per_out[sp] = {}
            if rising is not None:
                per_out[sp]["rising"] = rising
            if falling is not None:
                per_out[sp]["falling"] = falling
        if per_out:
            parsed["per_species"] = per_out
        out["outlet_sopdt"] = parsed
    return out


def parse_cape_request(data: dict[str, Any]) -> dict[str, Any]:
    """Validate and normalize a ``cape_request_v1`` object."""
    if not isinstance(data, dict):
        raise CapeContractError("request must be a JSON object")
    schema = str(data.get("schema") or "")
    if schema != REQUEST_SCHEMA:
        raise CapeContractError(f"schema must be {REQUEST_SCHEMA!r}, got {schema!r}")
    factory_id = validate_factory_id(data.get("factory_id"))

    config_path = data.get("config_path")
    config = data.get("config")
    if not config_path and not isinstance(config, dict):
        raise CapeContractError("provide config_path and/or inline config")

    feeds_raw = data.get("feeds")
    if not isinstance(feeds_raw, dict):
        raise CapeContractError("feeds must be an object")
    feeds: dict[str, Any] = {}
    for fid in FEED_IDS:
        parsed = _parse_feed(
            feeds_raw.get(fid), fid, required=fid in REQUIRED_FEEDS
        )
        if parsed is not None:
            feeds[fid] = parsed

    c_t = str(data.get("C_t_export") or "both").strip().lower()
    if c_t not in C_T_EXPORT_MODES:
        raise CapeContractError(
            f"C_t_export must be one of {C_T_EXPORT_MODES}, got {c_t!r}"
        )

    params = parse_unit_parameters(data.get("parameters"))

    return {
        "schema": REQUEST_SCHEMA,
        "factory_id": factory_id,
        "config_path": str(config_path) if config_path else None,
        "config": config if isinstance(config, dict) else None,
        "feeds": feeds,
        "parameters": params,
        "C_t_export": c_t,
        "balance_tol": float(data.get("balance_tol") or BALANCE_TOL_DEFAULT),
        "fail_if_unbalanced": bool(data.get("fail_if_unbalanced", False)),
    }


def load_cape_request(path: Path | str) -> dict[str, Any]:
    p = Path(path)
    data = json.loads(p.read_text(encoding="utf-8"))
    return parse_cape_request(data)


def dump_cape_json(payload: dict[str, Any], path: Path | str) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return p


def _format_balance_report(summary: dict[str, Any]) -> str:
    """DFM atom-ledger text (no dolfinx). Fermentation schemas are rejected."""
    schema = str(summary.get("schema") or "")
    if schema.startswith("fermentation"):
        raise CapeContractError(
            "balance_summary schema is fermentation; this PMC is DFM-only"
        )
    if "F_in_mol" not in summary or "atom_residual_mol" not in summary:
        raise CapeContractError("balance_summary is not a DFM atom ledger")
    lines = ["--- Mass / atom balance summary ---"]
    fin = summary["F_in_mol"]
    fout = summary["F_out_mol"]
    lines.append("  F_in  [mol]: " + ", ".join(f"{sp}={fin[sp]:.6g}" for sp in fin))
    lines.append("  F_out [mol]: " + ", ".join(f"{sp}={fout[sp]:.6g}" for sp in fout))
    h0 = summary.get("holdup_initial_mol") or {}
    hf = summary.get("holdup_final_mol") or {}
    lines.append(
        f"  N_CO2_ads: {h0.get('CO2_ads', 0):.6g} → {hf.get('CO2_ads', 0):.6g} mol"
    )
    ar = summary["atom_residual_mol"]
    arr = summary["atom_residual_relative"]
    lines.append(
        "  Atom residual R=ΔN−(Fin−Fout) [mol] (rel): "
        + ", ".join(f"{el}={ar[el]:.3e} ({arr[el]:.3e})" for el in ("C", "H", "O"))
    )
    conv = summary.get("conversion") or {}
    lines.append(
        f"  n_CH4={conv.get('n_CH4_produced_mol')} mol, "
        f"H2O/CH4={conv.get('H2O_over_CH4')}, "
        f"Y_CH4={conv.get('Y_CH4_per_CO2_in')}, "
        f"X_H2→CH4={conv.get('X_H2_to_CH4')}, "
        f"X_CO2_converted={conv.get('X_CO2_converted')}"
    )
    return "\n".join(lines)


def _horizon_s(run: dict[str, Any]) -> float:
    times = run.get("times_s") or []
    if not times:
        raise CapeContractError("run JSON has no times_s")
    t_end = float(times[-1])
    if t_end <= 0.0:
        raise CapeContractError("run horizon t_end must be > 0")
    return t_end


def cycle_average_product(
    F_out_mol: dict[str, float],
    t_horizon_s: float,
) -> dict[str, Any]:
    """Raw product: F_i^ss = F_i^out / t_horizon. SOPDT is not applied."""
    if t_horizon_s <= 0.0:
        raise CapeContractError("t_horizon_s must be > 0")
    f_ss = {sp: float(F_out_mol.get(sp, 0.0)) / t_horizon_s for sp in SPECIES}
    total = sum(max(v, 0.0) for v in f_ss.values())
    if total <= 0.0:
        y = {sp: 0.0 for sp in SPECIES}
    else:
        y = {sp: max(f_ss[sp], 0.0) / total for sp in SPECIES}
    return {
        "F_ss_mol_s": f_ss,
        "y": y,
        "y_cas": {SPECIES_TO_CAS[sp]: y[sp] for sp in SPECIES},
        "note": "raw cycle-average of bed outlet; SOPDT not applied",
    }


def build_cape_result(
    run: dict[str, Any],
    request: dict[str, Any] | None = None,
    *,
    T_K: float | None = None,
    P_Pa: float | None = None,
) -> dict[str, Any]:
    """
    Build ``cape_result_v1`` from a DFM ``fenics_run`` export.

    Product uses raw ``F_out_mol / t_end``. SOPDT series stay ``None`` here;
    ``to_cape_result`` (step 2) fills the report series from the config filter.
    Raw ``C(t)`` is copied when requested.
    """
    req = request or {
        "factory_id": ALLOWED_FACTORY_ID,
        "C_t_export": "both",
        "balance_tol": BALANCE_TOL_DEFAULT,
        "feeds": {},
    }
    if req.get("schema") == REQUEST_SCHEMA or "feeds" in req:
        if "factory_id" in req:
            validate_factory_id(req["factory_id"])

    meta = run.get("meta") or {}
    run_fid = str(meta.get("factory_id") or req.get("factory_id") or "")
    if run_fid and run_fid != ALLOWED_FACTORY_ID:
        validate_factory_id(run_fid)

    summary = run.get("balance_summary")
    from adr.cape.balance import (
        dfm_outputs_from_summary,
        enforce_balance_gate,
        require_dfm_ledger,
    )

    summary = require_dfm_ledger(summary)

    t_h = _horizon_s(run)
    fout = run.get("F_out_mol") or summary.get("F_out_mol") or {}
    product = cycle_average_product(fout, t_h)

    feeds = (req.get("feeds") or {}) if isinstance(req, dict) else {}
    ads = feeds.get("ads") or {}
    p_bed = float(P_Pa if P_Pa is not None else ads.get("P_Pa") or 101325.0)
    t_bed = float(T_K if T_K is not None else 593.15)
    product["T_K"] = t_bed
    product["P_Pa"] = p_bed

    tol = float(req.get("balance_tol") or BALANCE_TOL_DEFAULT)
    outputs = dfm_outputs_from_summary(summary, tol=tol)
    enforce_balance_gate(
        outputs, fail_if_unbalanced=bool(req.get("fail_if_unbalanced"))
    )
    ok = bool(outputs["balance_ok"])
    c_t = str(req.get("C_t_export") or "both")

    result: dict[str, Any] = {
        "schema": RESULT_SCHEMA,
        "factory_id": ALLOWED_FACTORY_ID,
        "t_horizon_s": t_h,
        "product": product,
        "outputs": outputs,
        "balance_ok": ok,
        "balance_tol": tol,
        "balance_summary": summary,
        "report_text": _format_balance_report(summary),
        "C_t_export": c_t,
    }

    times = run.get("times_s") or []
    raw_c = run.get("outlet_C_mol_m3") or {}
    if c_t in ("raw", "both") and times and raw_c:
        result["outlet_C_raw"] = {
            "times_s": [float(t) for t in times],
            "C_mol_m3": {
                sp: [float(v) for v in raw_c[sp]] for sp in SPECIES if sp in raw_c
            },
        }
    if c_t in ("sopdt", "both"):
        result["outlet_C_sopdt"] = None
        result["outlet_C_sopdt_note"] = (
            "SOPDT C(t) is a report series; filled in step 2 from outlet_sopdt_filter. "
            "Not used for product F_ss."
        )
    profiles = run.get("axial_profiles")
    if profiles:
        result["axial_profiles"] = profiles
    from adr.cape.reports import attach_reports

    attach_reports(result)
    return result
