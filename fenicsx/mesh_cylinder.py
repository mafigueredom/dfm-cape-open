"""Build axisymmetric (r, z) mesh: r in [0, R], z in [0, L]."""

from __future__ import annotations

import numpy as np
from mpi4py import MPI
from dolfinx import mesh


def create_axisymmetric_mesh(
    R_m: float,
    L_m: float,
    n_cells_r: int,
    n_cells_z: int,
    comm=MPI.COMM_WORLD,
):
    if min(R_m, L_m) <= 0:
        raise ValueError("R and L must be positive")
    if n_cells_r < 1 or n_cells_z < 1:
        raise ValueError("Need at least one cell in r and z")
    # dolfinx coordinates: x[0]=r, x[1]=z
    return mesh.create_rectangle(
        comm,
        [np.array([0.0, 0.0]), np.array([float(R_m), float(L_m)])],
        [int(n_cells_r), int(n_cells_z)],
        cell_type=mesh.CellType.triangle,
    )
