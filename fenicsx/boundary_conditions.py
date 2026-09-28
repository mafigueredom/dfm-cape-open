"""
Boundary facets and model-sheet BC reference (axisymmetric domain 0 ≤ r ≤ R, 0 ≤ z ≤ L).

Slide «Consideraciones en las Ecuaciones» — Fronteras / Condiciones iniciales
------------------------------------------------------------------------------

Momentum (Darcy):
  Entrada:  -n·(ρ u) = ρ U₀          [v1.1: prescribed u_z = Q/(ε π R²)]
  Salida:   P = P_op                  [implemented: update_pressure_profile → P(L) = P_op]
  Pared:    -n·(ρ u) = 0              [v1: no explicit wall constraint on u; axis r=0 is symmetric]

Mass (diluted species C_i):
  Entrada (COMSOL FluxDanckwerts, v1.2):
            n·(-D_z ∇C_i + u C_i) = n·(u C_{i,0})
            → weak form adds ∫_inlet r u_z (C - C_{i,0}) v ds  (u_z interstitial, +z)
  Entrada (v1.0 Dirichlet fallback):
            C_i = C_{i,0} m̂_i(t)
  Salida:   n·(D_z ∇C_i) = 0          [natural Neumann in diluted_species weak form]
  Pared:    n·(D_z ∇C_i) = 0          [natural on r = R; axis r = 0 via axisymmetric measure]

Initial:
  P = P_op;  C_N2 = P_op/(R T_op);  C_other = 0;  θ_CO2 = 0
"""

from __future__ import annotations

import numpy as np
from dolfinx import mesh

from fenicsx.properties import inlet_bc_mode

__all__ = [
    "inlet_bc_mode",
    "create_inlet_facet_tags",
    "locate_inlet_facets",
    "locate_outlet_facets",
    "locate_wall_facets",
]


def create_inlet_facet_tags(domain, inlet_facets, *, tag: int = 1):
    """Mesh tags for inlet facets (ds integrals)."""
    fdim = domain.topology.dim - 1
    values = np.full(len(inlet_facets), int(tag), dtype=np.int32)
    return mesh.meshtags(domain, fdim, inlet_facets, values)


def locate_inlet_facets(domain, fdim: int | None = None):
    """Facets at z = 0 (inlet)."""
    fdim = fdim if fdim is not None else domain.topology.dim - 1

    def on_inlet(x):
        return np.isclose(x[1], 0.0)

    return mesh.locate_entities_boundary(domain, fdim, on_inlet)


def locate_outlet_facets(domain, L_m: float, fdim: int | None = None):
    """Facets at z = L (outlet)."""
    fdim = fdim if fdim is not None else domain.topology.dim - 1

    def on_outlet(x):
        return np.isclose(x[1], float(L_m))

    return mesh.locate_entities_boundary(domain, fdim, on_outlet)


def locate_wall_facets(domain, R_m: float, fdim: int | None = None):
    """Facets at r = R (cylinder wall). Zero diffusive flux is natural (no extra BC object)."""
    fdim = fdim if fdim is not None else domain.topology.dim - 1

    def on_wall(x):
        return np.isclose(x[0], float(R_m))

    return mesh.locate_entities_boundary(domain, fdim, on_wall)
