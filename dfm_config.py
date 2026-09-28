"""Load full DFM configuration from JSON (schema 1.4)."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from geometry_config import GeometryConfig, MeshSpec, _mesh_from_payload, compute_mesh


@dataclass(frozen=True)
class DfmConfig:
    geometry: GeometryConfig
    gas_phase: dict[str, Any]
    solid_kinetics: dict[str, Any]
    energy_thermal: dict[str, Any]
    raw: dict[str, Any]


def dfm_config_from_dict(data: dict[str, Any]) -> DfmConfig:
    """Build ``DfmConfig`` from an in-memory JSON object (does not write disk)."""
    g = data["geometry"]
    b = data["bed_structure"]
    d_p = b.get("d_p_m")
    mesh = _mesh_from_payload(data["mesh"])
    geom = GeometryConfig(
        L_m=float(g["L_m"]),
        R_m=float(g["R_m"]),
        geometry_kind=str(g["kind"]),
        epsilon=float(b["epsilon"]),
        rho_s_kg_m3=float(b["rho_s_kg_m3"]),
        d_p_m=float(d_p) if d_p is not None else None,
        use_d_p_in_transport_correlations=bool(
            data.get("gas_phase", {})
            .get("diffusion", {})
            .get("use_wakao_smith", False)
        ),
        mesh=mesh,
        raw=data,
    )
    return DfmConfig(
        geometry=geom,
        gas_phase=dict(data["gas_phase"]),
        solid_kinetics=dict(data.get("solid_kinetics", {})),
        energy_thermal=dict(data.get("energy_thermal", {})),
        raw=data,
    )


def load_dfm_config(path: Path | str) -> DfmConfig:
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    return dfm_config_from_dict(data)


__all__ = [
    "DfmConfig",
    "dfm_config_from_dict",
    "load_dfm_config",
    "compute_mesh",
    "MeshSpec",
]
