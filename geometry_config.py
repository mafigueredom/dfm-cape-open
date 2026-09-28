"""Load and validate block-A geometry configuration (SI)."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class MeshSpec:
    Nz: int
    Nr: int
    dz_m: float
    dr_m: float


@dataclass(frozen=True)
class GeometryConfig:
    L_m: float
    R_m: float
    geometry_kind: str
    epsilon: float
    rho_s_kg_m3: float
    d_p_m: float | None
    use_d_p_in_transport_correlations: bool
    mesh: MeshSpec
    raw: dict[str, Any]

    @property
    def D_m(self) -> float:
        return 2.0 * self.R_m


def _mesh_from_payload(mesh: dict[str, Any]) -> MeshSpec:
    nodes = mesh["nodes"]
    eff = mesh["spacing_effective_m"]
    return MeshSpec(
        Nz=int(nodes["Nz"]),
        Nr=int(nodes["Nr"]),
        dz_m=float(eff["dz"]),
        dr_m=float(eff["dr"]),
    )


def load_geometry_config(path: Path | str) -> GeometryConfig:
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    g = data["geometry"]
    b = data["bed_structure"]
    d_p = b.get("d_p_m")
    return GeometryConfig(
        L_m=float(g["L_m"]),
        R_m=float(g["R_m"]),
        geometry_kind=str(g["kind"]),
        epsilon=float(b["epsilon"]),
        rho_s_kg_m3=float(b["rho_s_kg_m3"]),
        d_p_m=float(d_p) if d_p is not None else None,
        use_d_p_in_transport_correlations=bool(b.get("use_d_p_in_transport_correlations", False)),
        mesh=_mesh_from_payload(data["mesh"]),
        raw=data,
    )


def compute_mesh(L_m: float, R_m: float, dz_m: float, dr_m: float) -> MeshSpec:
    """Same rule as web UI: uniform grid, boundaries at z=0,L and r=0,R."""
    if min(L_m, R_m, dz_m, dr_m) <= 0:
        raise ValueError("L, R, dz, dr must be positive")
    Nz = max(2, int(L_m // dz_m) + 1)
    Nr = max(2, int(R_m // dr_m) + 1)
    return MeshSpec(
        Nz=Nz,
        Nr=Nr,
        dz_m=L_m / (Nz - 1),
        dr_m=R_m / (Nr - 1),
    )
