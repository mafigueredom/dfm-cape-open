"""
Darcy flow — mass continuity + Darcy velocity law.

Mass continuity (porosity ε_cat):
    ε_cat ∂ρ_f/∂t + ∇·(ρ_f u) = Σ R_i

Darcy law (model sheet):
    u = (𝒦 / μ_f) ∇P

Pressure profile: v1 uses axial linear P(z) consistent with cycle interstitial v_z until
a full pressure Poisson equation is added in a later sheet.
"""

from __future__ import annotations

import numpy as np
import ufl
from dolfinx import fem, default_scalar_type
from dolfinx.fem.petsc import LinearProblem

from dfm_config import DfmConfig
from fenicsx.properties import (
    darcy_permeability_m2,
    fluid_viscosity_Pa_s,
    molar_mass_kg,
    species_list,
    total_pressure_Pa,
)


def axisymmetric_divergence(u: ufl.core.expr.Expr, r: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
    """∇·u in (r, z) axisymmetric reduction: (1/r)∂(r u_r)/∂r + ∂u_z/∂z."""
    u_r = u[0]
    u_z = u[1]
    return (1.0 / (r + 1e-30)) * ufl.Dx(r * u_r, 0) + ufl.Dx(u_z, 1)


def mass_flux_divergence(rho: ufl.core.expr.Expr, u: ufl.core.expr.Expr, r: ufl.core.expr.Expr) -> ufl.core.expr.Expr:
    """∇·(ρ u) for the continuity equation."""
    return axisymmetric_divergence(rho * u, r)


def darcy_velocity(P: fem.Function, cfg: DfmConfig) -> ufl.core.expr.Expr:
    """u = (𝒦/μ_f) ∇P as on the model sheet (sign convention without leading minus)."""
    K = darcy_permeability_m2(cfg)
    mu = fluid_viscosity_Pa_s(cfg)
    return (K / mu) * ufl.grad(P)


def update_pressure_profile(P_fn: fem.Function, v_z: float, cfg: DfmConfig) -> None:
    """
    Axial pressure with outlet BC P = P_op (model sheet — Salida: P = P_op).

    Darcy relation u ≈ (𝒦/μ_f) ∂P/∂z  →  P(z) = P_op + (μ_f/𝒦) v_z (z − L).
    """
    K = darcy_permeability_m2(cfg)
    mu = fluid_viscosity_Pa_s(cfg)
    P_op = total_pressure_Pa(cfg)
    L_m = float(cfg.geometry.L_m)
    coords = P_fn.function_space.tabulate_dof_coordinates()
    z = coords[:, 1]
    dP_dz = (mu / K) * float(v_z)
    P_fn.x.array[:] = P_op + dP_dz * (z - L_m)


def update_prescribed_velocity(u_fn: fem.Function, v_z_m_s: float) -> None:
    """
    Uniform axial interstitial velocity (COMSOL POW_V07: u = dl.w/e_dfm, u_r = 0).

    v_z_m_s = Q / (epsilon * pi * R^2).

    Vector Lagrange dofs are interleaved per node: [u_r0, u_z0, u_r1, u_z1, ...].
    """
    arr = u_fn.x.array
    arr[0::2] = 0.0
    arr[1::2] = float(v_z_m_s)


def velocity_model(cfg: DfmConfig) -> str:
    """prescribed_u0 | darcy_pressure (default)."""
    vel = cfg.gas_phase.get("velocity") or {}
    mom = cfg.gas_phase.get("momentum") or {}
    raw = vel.get("model") or mom.get("velocity_model") or "darcy_pressure"
    return str(raw).strip().lower()


def update_velocity_field(
    u_fn: fem.Function,
    P_fn: fem.Function,
    v_z_m_s: float,
    cfg: DfmConfig,
) -> None:
    if velocity_model(cfg) in ("prescribed_u0", "prescribed", "u0", "comsol_vel_d"):
        update_prescribed_velocity(u_fn, v_z_m_s)
    else:
        update_velocity_from_darcy(u_fn, P_fn, cfg)


def update_velocity_from_darcy(u_fn: fem.Function, P_fn: fem.Function, cfg: DfmConfig) -> None:
    """L² projection of u = (𝒦/μ_f) ∇P into the vector Lagrange space."""
    V_u = u_fn.function_space
    u = ufl.TrialFunction(V_u)
    v = ufl.TestFunction(V_u)
    u_expr = darcy_velocity(P_fn, cfg)
    a = ufl.inner(u, v) * ufl.dx
    L = ufl.inner(u_expr, v) * ufl.dx
    problem = LinearProblem(
        a,
        L,
        bcs=[],
        petsc_options_prefix="darcy_u_",
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    )
    uh = problem.solve()
    u_fn.x.array[:] = uh.x.array[:]


def update_fluid_density(rho_fn: fem.Function, C_fns: dict[str, fem.Function], cfg: DfmConfig) -> None:
    """Pointwise ρ_f = Σ_i C_i M_i from current gas concentrations (vectorized)."""
    rho = rho_fn.x.array
    rho[:] = 0.0
    for sp in species_list(cfg):
        rho += molar_mass_kg(sp, cfg) * C_fns[sp].x.array


def continuity_weak_form(
    rho: ufl.TrialFunction,
    rho_n: fem.Function,
    u: fem.Function,
    source: fem.Function,
    *,
    epsilon: float,
    dt: float,
    r_coord: ufl.core.expr.Expr,
    test: ufl.TestFunction,
):
    """
    Backward-Euler weak form of ε ∂ρ/∂t + ∇·(ρ u) = S with integration-by-parts on advection.

    Find ρ^{n+1} such that for all v:
        ∫ r ε (ρ−ρ^n)/dt v dx  − ∫ r (ρ^n u)·∇v dx  = ∫ r S v dx
    """
    eps_c = fem.Constant(rho_n.function_space.mesh, default_scalar_type(epsilon))
    dt_c = fem.Constant(rho_n.function_space.mesh, default_scalar_type(dt))
    r = r_coord
    div_term = ufl.dot(rho_n * u, ufl.grad(test)) * r * ufl.dx
    a = eps_c / dt_c * rho * test * r * ufl.dx
    L = eps_c / dt_c * rho_n * test * r * ufl.dx + source * test * r * ufl.dx - div_term
    return a, L


class ContinuitySolver:
    """
    Cached implicit continuity step (same weak form as solve_continuity_implicit).

    LHS (ε/dt mass matrix) depends only on dt — LU is reused until dt changes.
    RHS (ρ_n, u, source) is reassembled every step.
    """

    def __init__(
        self,
        rho_n: fem.Function,
        rho_new: fem.Function,
        u: fem.Function,
        source: fem.Function,
        *,
        epsilon: float,
        r_coord: ufl.core.expr.Expr,
    ) -> None:
        from fenicsx.diluted_species import CachedLinearSolver

        V = rho_n.function_space
        mesh = V.mesh
        rho = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)
        r = r_coord
        self._rho_new = rho_new
        self._dt = fem.Constant(mesh, default_scalar_type(1.0))
        self._dt_cached: float | None = None
        eps_c = fem.Constant(mesh, default_scalar_type(epsilon))
        div_term = ufl.dot(rho_n * u, ufl.grad(v)) * r * ufl.dx
        a = eps_c / self._dt * rho * v * r * ufl.dx
        L = eps_c / self._dt * rho_n * v * r * ufl.dx + source * v * r * ufl.dx - div_term
        self._uh = fem.Function(V)
        self._solver = CachedLinearSolver(
            a, L, bcs=[], petsc_prefix="darcy_rho_", u_out=self._uh
        )

    def mark_matrix_dirty(self) -> None:
        self._solver.mark_matrix_dirty()

    def step(self, dt: float) -> None:
        if self._dt_cached is None or abs(float(dt) - self._dt_cached) > 1e-15:
            self._dt.value = dt
            self._dt_cached = float(dt)
            self._solver.mark_matrix_dirty()
        self._solver.solve()
        self._rho_new.x.array[:] = self._uh.x.array[:]


def solve_continuity_implicit(
    rho_n: fem.Function,
    rho_new: fem.Function,
    u: fem.Function,
    source: fem.Function,
    *,
    epsilon: float,
    dt: float,
    r_coord: ufl.core.expr.Expr,
    bcs: list | None = None,
) -> fem.Function:
    """Implicit Galerkin step for ρ_f (linearized with ρ^n in advection flux)."""
    V = rho_n.function_space
    rho = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a, L = continuity_weak_form(
        rho, rho_n, u, source, epsilon=epsilon, dt=dt, r_coord=r_coord, test=v
    )
    problem = LinearProblem(
        a,
        L,
        bcs=bcs or [],
        petsc_options_prefix="darcy_rho_",
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    )
    uh = problem.solve()
    rho_new.x.array[:] = uh.x.array[:]
    return rho_new
