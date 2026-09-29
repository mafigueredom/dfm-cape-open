"""Serialize / load FEniCS reactor runs for offline comparison (v1.4+)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from adr.model.state import FieldStateHistory as ReactorState


def reactor_state_to_dict(
    state: ReactorState,
    *,
    species: list[str],
    config_path: str | Path | None = None,
    meta: dict[str, Any] | None = None,
) -> dict[str, Any]:
    out: dict[str, Any] = {
        "schema": "fenics_run_v1.1",
        "species": list(species),
        "times_s": [float(t) for t in state.times],
        "outlet_C_mol_m3": {sp: [float(v) for v in state.outlet_means[sp]] for sp in species},
        "inlet_C_mol_m3": {
            sp: [float(v) for v in state.inlet_effective.get(sp, [])] for sp in species
        },
        "phase_id": list(state.phase_id),
    }
    if state.Q_m3_s:
        out["Q_m3_s"] = [float(q) for q in state.Q_m3_s]
    if state.holdup_mol:
        out["holdup_mol"] = {
            k: [float(v) for v in series] for k, series in state.holdup_mol.items()
        }
    if state.F_in_mol:
        out["F_in_mol"] = {sp: float(state.F_in_mol.get(sp, 0.0)) for sp in species}
    if state.F_out_mol:
        out["F_out_mol"] = {sp: float(state.F_out_mol.get(sp, 0.0)) for sp in species}
    if state.balance_summary is not None:
        out["balance_summary"] = state.balance_summary
    if state.axial_profiles:
        out["axial_profiles"] = state.axial_profiles
    if config_path is not None:
        out["config"] = str(config_path)
    if meta:
        out["meta"] = meta
    return out


def save_reactor_state_json(
    path: str | Path,
    state: ReactorState,
    *,
    species: list[str],
    config_path: str | Path | None = None,
    meta: dict[str, Any] | None = None,
) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    payload = reactor_state_to_dict(
        state, species=species, config_path=config_path, meta=meta
    )
    p.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return p


def load_fenics_run_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))
