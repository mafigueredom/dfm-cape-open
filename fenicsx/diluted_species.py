"""
Diluted species transport (model sheet — Masa: Transp. Especies Diluidas).

Strong form per species i:
    ∂C_i/∂t − D_z ∇²C_i + u·∇C_i = R_i

Axial dispersion:
    D_z = 0.7 D_m + 0.5 d_p |u|

Boundary conditions (weak / natural; see `fenicsx/boundary_conditions.py`):

    Entrada (v1.2 Danckwerts / COMSOL FluxDanckwerts):
        n·(-D_z ∇C_i + u C_i) = n·(u C_{i,0})
        → adds ∫_Γ_in r u_z (C − C_{i,0}) v ds to the weak form
    Entrada (v1.0 Dirichlet fallback):
        C_i = C_{i,0} m̂_i(t)
    Salida (masa):   n·(D_z ∇C_i) = 0   [natural in weak form]
    Pared (masa):    n·(D_z ∇C_i) = 0   [natural on r = R]
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import ufl
from dolfinx import fem, default_scalar_type
from dolfinx.fem import petsc as fem_petsc
from petsc4py import PETSc


@dataclass(frozen=True)
class DanckwertsInletBC:
    """Inlet Danckwerts BC data for one species step."""

    facet_tags: object
    subdomain_id: int
    u_z_m_s: float
    C_in: fem.Function


class CachedLinearSolver:
    """
    Assemble-once linear solver: LU of A is reused until mark_matrix_dirty().

    Each step always reassembles the RHS vector; the left-hand-side matrix is
    reassembled and refactorized only when coefficients that enter A change
    (dt, velocity, dispersion, …).
    """

    def __init__(
        self,
        a_form,
        L_form,
        *,
        bcs: list,
        petsc_prefix: str,
        u_out: fem.Function,
    ) -> None:
        self._a = fem.form(a_form)
        self._L = fem.form(L_form)
        self._bcs = bcs
        self._u = u_out
        V = u_out.function_space
        self._A = fem_petsc.create_matrix(self._a)
        self._b = fem_petsc.create_vector(V)
        self._ksp = PETSc.KSP().create(self._A.getComm())
        self._ksp.setOperators(self._A)
        self._ksp.setType("preonly")
        pc = self._ksp.getPC()
        pc.setType("lu")
        self._ksp.setOptionsPrefix(petsc_prefix)
        self._ksp.setFromOptions()
        self._matrix_dirty = True

    def mark_matrix_dirty(self) -> None:
        self._matrix_dirty = True

    def solve(self) -> None:
        if self._matrix_dirty:
            self._A.zeroEntries()
            fem_petsc.assemble_matrix(self._A, self._a, bcs=self._bcs)
            self._A.assemble()
            self._ksp.setOperators(self._A)
            self._matrix_dirty = False

        self._b.set(0.0)
        fem_petsc.assemble_vector(self._b, self._L)
        fem_petsc.apply_lifting(self._b, [self._a], bcs=[self._bcs])
        self._b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        fem_petsc.set_bc(self._b, self._bcs)
        self._ksp.solve(self._b, self._u.x.petsc_vec)
        self._u.x.scatter_forward()

    def destroy(self) -> None:
        self._ksp.destroy()
        self._A.destroy()
        self._b.destroy()


class SpeciesTransportSolver:
    """
    Cached backward-Euler transport solver.

    Forms/PETSc objects are built once. The LHS matrix (and its LU) is reused
    across steps until dt or other A-coefficients change; call
    mark_matrix_dirty() after updating D_z / u, or pass a new dt via step().
    """

    def __init__(
        self,
        *,
        C_n: fem.Function,
        D_z: fem.Function,
        u: fem.Function,
        source: fem.Function,
        r_coord: ufl.core.expr.Expr,
        bcs: list,
        petsc_prefix: str,
        danckwerts: DanckwertsInletBC | None = None,
        supg: bool = False,
    ) -> None:
        domain = C_n.function_space.mesh
        V = C_n.function_space
        C = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)
        r = r_coord
        self._C_n = C_n
        self._uh = fem.Function(V)
        self._dt = fem.Constant(domain, default_scalar_type(1.0))
        self._dt_cached: float | None = None

        D_eff = D_z
        if supg:
            h = ufl.CellDiameter(domain)
            u_mag = ufl.sqrt(ufl.dot(u, u) + 1e-30)
            D_eff = D_z + 0.5 * u_mag * h

        a = (
            r * C * v * ufl.dx
            + self._dt * r * D_eff * ufl.dot(ufl.grad(C), ufl.grad(v)) * ufl.dx
            + self._dt * r * ufl.dot(u, ufl.grad(C)) * v * ufl.dx
        )
        L = r * C_n * v * ufl.dx + self._dt * r * source * v * ufl.dx

        self._u_z: fem.Constant | None = None
        if danckwerts is not None:
            ds_in = ufl.Measure(
                "ds",
                domain=domain,
                subdomain_data=danckwerts.facet_tags,
                subdomain_id=int(danckwerts.subdomain_id),
            )
            self._u_z = fem.Constant(domain, default_scalar_type(danckwerts.u_z_m_s))
            a += self._dt * r * self._u_z * C * v * ds_in
            L += self._dt * r * self._u_z * danckwerts.C_in * v * ds_in

        self._solver = CachedLinearSolver(
            a, L, bcs=bcs, petsc_prefix=petsc_prefix, u_out=self._uh
        )

    def mark_matrix_dirty(self) -> None:
        self._solver.mark_matrix_dirty()

    def step(self, dt: float, u_z_m_s: float | None = None) -> None:
        """Advance C_n by one backward-Euler step of size dt (clips negatives)."""
        if self._dt_cached is None or abs(float(dt) - self._dt_cached) > 1e-15:
            self._dt.value = dt
            self._dt_cached = float(dt)
            self._solver.mark_matrix_dirty()
        if self._u_z is not None and u_z_m_s is not None:
            if abs(float(self._u_z.value) - float(u_z_m_s)) > 1e-15:
                self._u_z.value = u_z_m_s
                self._solver.mark_matrix_dirty()
        self._solver.solve()
        self._C_n.x.array[:] = np.maximum(self._uh.x.array[:], 0.0)


def update_D_z_field(
    D_z_fn: fem.Function,
    u_fn: fem.Function,
    D_m: float,
    d_p: float,
    *,
    alpha: float = 0.7,
    beta: float = 0.5,
) -> None:
    """Pointwise D_z = alpha D_m + beta d_p |u|."""
    u_arr = u_fn.x.array
    n_nodes = len(D_z_fn.x.array)
    u_mag = np.zeros(n_nodes, dtype=default_scalar_type)
    n_pairs = len(u_arr) // 2
    if n_pairs == n_nodes:
        u_mag = np.sqrt(u_arr[0::2] ** 2 + u_arr[1::2] ** 2)
    else:
        u_mag[:] = np.abs(u_arr[:n_nodes])
    D_z_fn.x.array[:] = alpha * float(D_m) + beta * float(d_p) * u_mag


def solve_diluted_species_step(
    *,
    C_n: fem.Function,
    D_z: fem.Function,
    u: fem.Function,
    source: fem.Function,
    dt: float,
    r_coord: ufl.core.expr.Expr,
    bcs: list,
    petsc_prefix: str,
    danckwerts: DanckwertsInletBC | None = None,
    supg: bool = False,
    sink: fem.Function | None = None,
) -> None:
    """
    One-shot backward-Euler step (forms rebuilt each call). Prefer
    SpeciesTransportSolver in the time loop.
    """
    domain = C_n.function_space.mesh
    V = C_n.function_space
    C = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    r = r_coord
    dt_c = fem.Constant(domain, default_scalar_type(dt))

    D_eff = D_z
    if supg:
        h = ufl.CellDiameter(domain)
        u_mag = ufl.sqrt(ufl.dot(u, u) + 1e-30)
        D_eff = D_z + 0.5 * u_mag * h

    a = (
        r * C * v * ufl.dx
        + dt_c * r * D_eff * ufl.dot(ufl.grad(C), ufl.grad(v)) * ufl.dx
        + dt_c * r * ufl.dot(u, ufl.grad(C)) * v * ufl.dx
    )
    if sink is not None:
        a += dt_c * r * sink * C * v * ufl.dx
    L = r * C_n * v * ufl.dx + dt_c * r * source * v * ufl.dx

    if danckwerts is not None:
        ds_in = ufl.Measure(
            "ds",
            domain=domain,
            subdomain_data=danckwerts.facet_tags,
            subdomain_id=int(danckwerts.subdomain_id),
        )
        u_z = fem.Constant(domain, default_scalar_type(danckwerts.u_z_m_s))
        a += dt_c * r * u_z * C * v * ds_in
        L += dt_c * r * u_z * danckwerts.C_in * v * ds_in

    uh = fem.Function(V)
    solver = CachedLinearSolver(a, L, bcs=bcs, petsc_prefix=petsc_prefix, u_out=uh)
    solver.solve()
    C_n.x.array[:] = np.maximum(uh.x.array[:], 0.0)
    solver.destroy()
