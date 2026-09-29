"""
FEniCSx runtime backend for the methanation / DFM problem family.

Holds mesh, function spaces, cached transport/continuity solvers, and
step helpers invoked by ``AdrSolver``. Physics formulation lives on
``AbstractProblem`` terms; this class is the numerical substrate.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
from mpi4py import MPI
from dolfinx import fem, default_scalar_type

from dfm_config import DfmConfig
from kinetics import (
    adsorption_exchange_arrays,
    csc_qe_evaluator,
    k_ads_si,
    make_mass_source_evaluator,
    make_reaction_exchange_arrays,
)
from fenicsx.boundary_conditions import (
    create_inlet_facet_tags,
    locate_inlet_facets,
)
from fenicsx.mesh_cylinder import create_axisymmetric_mesh
from fenicsx.properties import (
    D_z_species,
    interstitial_velocity,
    q_e_const_mol_kg,
    rho_a_bed_kg_m3,
    rho_a_over_epsilon,
    superficial_velocity,
    total_pressure_Pa,
)
from fenicsx.darcy_mass_balance import (
    ContinuitySolver,
    update_fluid_density,
    update_pressure_profile,
    update_velocity_field,
)
from fenicsx.diluted_species import DanckwertsInletBC, SpeciesTransportSolver
from fenicsx.wilke_transport import D_m_mode, make_wilke_dm_evaluator
from fenicsx.balance_diagnostics import (
    HoldupAssembler,
    balance_enabled,
    balance_every_n_steps,
)
from adr.model.conditions import InletData
from adr.model.problem import AbstractProblem

if TYPE_CHECKING:
    from adr.factories.methanation_dfm import DfmBoundaryConditions


class DfmFenicsRuntime:
    """FEM state + local step operations for methanation_dfm."""

    def __init__(self, cfg: DfmConfig, problem: AbstractProblem):
        self.cfg = cfg
        self.problem = problem
        self.comm = MPI.COMM_WORLD
        self.species = list(problem.species.names)

        g = cfg.geometry
        m = cfg.geometry.mesh
        ncr = max(1, m.Nr - 1)
        ncz = max(1, m.Nz - 1)
        self.domain = create_axisymmetric_mesh(g.R_m, g.L_m, ncr, ncz, comm=self.comm)
        self.V = fem.functionspace(self.domain, ("Lagrange", 1))

        bcs = problem.bcs
        self._scheduler = problem.aux.scheduler
        self._Q_default = float(problem.aux.extras.get("Q_default") or 0.0)
        if self._Q_default <= 0:
            raise ValueError("gas_phase.velocity.Q_m3_s_at_T_ref required.")

        self.C_in = self._scheduler.inlet_at(0.0)
        self._meth_active = self._scheduler.methanation_active_at(0.0)
        self._ads_active = self._scheduler.adsorption_active_at(0.0)
        self._current_phase_id = self._scheduler.phase_at(0.0).id
        self.v_z = interstitial_velocity(
            cfg, self._scheduler.Q_m3_s_at(0.0, self._Q_default)
        )
        if g.d_p_m is None or float(g.d_p_m) <= 0:
            raise ValueError("bed_structure.d_p_m required.")
        self._d_p = float(g.d_p_m)

        ic = problem.ics.gas_concentrations()
        q_init = problem.ics.solid_loading()
        self.C_n = {sp: fem.Function(self.V) for sp in self.species}
        self.q_n = fem.Function(self.V)
        self.theta_n = fem.Function(self.V)
        self.q_n.x.array[:] = q_init
        self._qe_ref = q_e_const_mol_kg(cfg)
        self.theta_n.x.array[:] = q_init / self._qe_ref
        for sp in self.species:
            self.C_n[sp].x.array[:] = ic[sp]

        self._D_z = {sp: fem.Function(self.V) for sp in self.species}
        self._S_species = {sp: fem.Function(self.V) for sp in self.species}
        self._qe_fn = csc_qe_evaluator(cfg)
        self._k_ads = k_ads_si(cfg)
        self._pref = rho_a_over_epsilon(cfg)
        self._reaction_step_arrays = make_reaction_exchange_arrays(cfg)
        self._sum_Ri = make_mass_source_evaluator(cfg)
        self.zero_source = np.zeros(
            len(self.q_n.x.array), dtype=default_scalar_type
        )

        tr = cfg.gas_phase.get("transport") or {}
        self._disp_alpha = float(tr.get("dispersion_alpha", 0.7))
        self._disp_beta = float(tr.get("dispersion_beta", 0.5))
        self._dz_mode = str(tr.get("D_z_mode", "dispersion_formula")).strip()
        self._wilke_eval = None
        if self._dz_mode == "dispersion_formula" and D_m_mode(cfg) == "wilke":
            self._wilke_eval = make_wilke_dm_evaluator(cfg, self.species)

        coords = self.V.tabulate_dof_coordinates()
        mask = np.isclose(coords[:, 1], cfg.geometry.L_m, rtol=0, atol=1e-9)
        self._outlet_mask = mask if np.any(mask) else None
        if self._outlet_mask is not None:
            w = np.maximum(coords[mask, 0], 1e-12)
            self._outlet_w = w / np.sum(w)

        self._build_dirichlet_scaffold()
        self._inlet_bc_mode = bcs.inlet_bc_mode()
        self._inlet_facet_tags = None
        if self._inlet_bc_mode == "danckwerts":
            fdim = self.domain.topology.dim - 1
            inlet_facets = locate_inlet_facets(self.domain, fdim)
            self._inlet_facet_tags = create_inlet_facet_tags(self.domain, inlet_facets)

        import ufl

        self._r_coord = ufl.SpatialCoordinate(self.domain)[0]
        self._epsilon = float(cfg.geometry.epsilon)
        self.rho_f_n = fem.Function(self.V)
        self.rho_f_new = fem.Function(self.V)
        self._mass_source = fem.Function(self.V)
        self.P_n = fem.Function(self.V)
        self.P_n.x.array[:] = total_pressure_Pa(cfg)
        self.u = fem.functionspace(
            self.domain, ("Lagrange", 1, (self.domain.geometry.dim,))
        )
        self.u_n = fem.Function(self.u)
        update_fluid_density(self.rho_f_n, self.C_n, cfg)
        self._vel_d = superficial_velocity(cfg, self._Q_default)
        self._last_v_z: float | None = None
        self._C_snap = {sp: self.C_n[sp].x.array.copy() for sp in self.species}
        self.transport_solvers: dict[str, SpeciesTransportSolver] = {}
        self._continuity_solver: ContinuitySolver | None = None

        # Seed velocity / D_z / inlet before building cached solvers
        self.apply_inlet_schedule(0.0, 0.0, bcs)  # type: ignore[arg-type]

        for sp in self.species:
            dkw = None
            if self._inlet_bc_mode == "danckwerts":
                dkw = DanckwertsInletBC(
                    facet_tags=self._inlet_facet_tags,
                    subdomain_id=1,
                    u_z_m_s=self.v_z,
                    C_in=self._C_bc[sp],
                )
            self.transport_solvers[sp] = SpeciesTransportSolver(
                C_n=self.C_n[sp],
                D_z=self._D_z[sp],
                u=self.u_n,
                source=self._S_species[sp],
                r_coord=self._r_coord,
                bcs=self._bcs[sp],
                petsc_prefix=f"diluted_{sp}_",
                danckwerts=dkw,
            )
        self._continuity_solver = ContinuitySolver(
            self.rho_f_n,
            self.rho_f_new,
            self.u_n,
            self._mass_source,
            epsilon=self._epsilon,
            r_coord=self._r_coord,
        )
        self.balance_on = balance_enabled(cfg)
        self.balance_every = balance_every_n_steps(cfg)
        self.holdup = (
            HoldupAssembler(
                self.domain,
                self._r_coord,
                epsilon=self._epsilon,
                rho_a=rho_a_bed_kg_m3(cfg),
                comm=self.comm,
            )
            if self.balance_on
            else None
        )

    def _build_dirichlet_scaffold(self) -> None:
        fdim = self.domain.topology.dim - 1
        inlet_facets = locate_inlet_facets(self.domain, fdim)
        self._inlet_dofs = fem.locate_dofs_topological(self.V, fdim, inlet_facets)
        self._C_bc = {sp: fem.Function(self.V) for sp in self.species}
        for sp in self.species:
            self._C_bc[sp].x.array[:] = self.C_in[sp]
        if self.problem.bcs.inlet_bc_mode() == "danckwerts":
            self._bcs = {sp: [] for sp in self.species}
        else:
            self._bcs = {
                sp: [fem.dirichletbc(self._C_bc[sp], self._inlet_dofs)]
                for sp in self.species
            }

    def _mark_transport_dirty(self) -> None:
        for s in self.transport_solvers.values():
            s.mark_matrix_dirty()
        if self._continuity_solver is not None:
            self._continuity_solver.mark_matrix_dirty()

    def apply_inlet_schedule(
        self, t: float, dt: float, bcs: DfmBoundaryConditions
    ) -> InletData:
        phase = bcs.scheduler.phase_at(t)
        if phase.id != bcs._last_phase_id:
            bcs._sopdt.begin_phase(t)
            bcs._last_phase_id = phase.id
        self._current_phase_id = phase.id
        self._meth_active = phase.methanation_active
        self._ads_active = phase.adsorption_active
        Q = bcs.scheduler.Q_m3_s_at(t, bcs._Q_default)
        v_z = interstitial_velocity(self.cfg, Q)
        vel_changed = self._last_v_z is None or abs(v_z - self._last_v_z) > 1e-15
        if vel_changed:
            self.v_z = v_z
            self._vel_d = superficial_velocity(self.cfg, Q)
            update_pressure_profile(self.P_n, self.v_z, self.cfg)
            update_velocity_field(self.u_n, self.P_n, self.v_z, self.cfg)
            self._last_v_z = v_z
            self._mark_transport_dirty()
        if self.update_D_z_fields(force=vel_changed) and self.transport_solvers:
            self._mark_transport_dirty()
        if dt > 0:
            cin_eff = bcs._sopdt.step(t, dt)
        else:
            cin_eff = bcs._sopdt.step(t, 1e-6)
        for sp in self.species:
            self._C_bc[sp].x.array[:] = cin_eff[sp]
        self.C_in = dict(cin_eff)
        return InletData(
            C_mol_m3=dict(cin_eff),
            Q_m3_s=Q,
            phase_id=phase.id,
            adsorption_active=bool(phase.adsorption_active),
            methanation_active=bool(phase.methanation_active),
        )

    def schedule_flags(self) -> tuple[float, float]:
        return float(self._ads_active), float(self._meth_active)

    def outlet_mean(self, C: fem.Function) -> float:
        if self._outlet_mask is None:
            return float(np.mean(C.x.array))
        return float(C.x.array[self._outlet_mask] @ self._outlet_w)

    def axial_profile(self) -> dict:
        """Radial-mean C_i(z) and q(z). q is solid CO2 loading [mol/kg]."""
        n_local = int(self.V.dofmap.index_map.size_local)
        coords = np.asarray(self.V.tabulate_dof_coordinates()[:n_local], dtype=float)
        r = coords[:, 0]
        z = coords[:, 1]
        q = np.asarray(self.q_n.x.array[:n_local], dtype=float)
        conc = {
            sp: np.asarray(self.C_n[sp].x.array[:n_local], dtype=float)
            for sp in self.species
        }
        if self.comm.size > 1:
            parts = self.comm.allgather((r, z, q, conc))
            r = np.concatenate([p[0] for p in parts])
            z = np.concatenate([p[1] for p in parts])
            q = np.concatenate([p[2] for p in parts])
            conc = {
                sp: np.concatenate([p[3][sp] for p in parts]) for sp in self.species
            }
        return _radial_mean_along_z(r, z, q, conc, self.species)

    def update_D_z_fields(self, *, force: bool = False) -> bool:
        vel_d = max(float(getattr(self, "_vel_d", self.v_z * self._epsilon)), 0.0)
        if self._wilke_eval is not None and not force:
            metric = 0.0
            for sp in self.species:
                arr = self.C_n[sp].x.array
                prev = self._C_snap[sp]
                scale = max(float(np.max(np.abs(arr))), 1e-9)
                metric = max(metric, float(np.max(np.abs(arr - prev))) / scale)
            if metric < 0.01:
                return False
        if self._wilke_eval is not None:
            Dm = self._wilke_eval({sp: self.C_n[sp].x.array for sp in self.species})
            for sp in self.species:
                self._D_z[sp].x.array[:] = (
                    self._disp_alpha * Dm[sp] + self._disp_beta * self._d_p * vel_d
                )
                self._C_snap[sp][:] = self.C_n[sp].x.array
            return True
        if force or self._last_v_z is None:
            for sp in self.species:
                self._D_z[sp].x.array[:] = D_z_species(
                    sp, self.v_z, self.cfg, C=None, vel_d_m_s=vel_d
                )
            return True
        return False

    def species_step(self, sp: str, dt: float, sources: np.ndarray) -> None:
        self._S_species[sp].x.array[:] = sources
        self.transport_solvers[sp].step(dt, u_z_m_s=self.v_z)

    def solid_gas_exchange(self, dt: float) -> None:
        s_co2, s_h2 = self.schedule_flags()
        q = self.q_n.x.array
        if s_co2 > 0 and self._k_ads > 0:
            C = self.C_n["CO2"].x.array
            Cn, qn = adsorption_exchange_arrays(
                C, q, dt, k_ads=self._k_ads, pref=self._pref, qe=self._qe_fn
            )
            C[:] = Cn
            q[:] = qn
        if s_h2 > 0:
            H = self.C_n["H2"].x.array
            Hn, qn, d_ch4, d_h2o, d_co2 = self._reaction_step_arrays(H, q, dt)
            H[:] = Hn
            q[:] = qn
            self.C_n["CH4"].x.array[:] += d_ch4
            self.C_n["H2O"].x.array[:] += d_h2o
            self.C_n["CO2"].x.array[:] += d_co2
        self.theta_n.x.array[:] = q / self._qe_ref

    def continuity_step(self, dt: float) -> None:
        s_co2, s_h2 = self.schedule_flags()
        c_co2 = self.C_n["CO2"].x.array if "CO2" in self.C_n else self.zero_source
        c_h2 = self.C_n["H2"].x.array if "H2" in self.C_n else self.zero_source
        self._mass_source.x.array[:] = self._sum_Ri(
            c_co2, c_h2, self.q_n.x.array, s_co2, s_h2
        )
        assert self._continuity_solver is not None
        self._continuity_solver.step(dt)
        self.rho_f_n.x.array[:] = self.rho_f_new.x.array[:]

    def update_density(self) -> None:
        update_fluid_density(self.rho_f_n, self.C_n, self.cfg)


def _radial_mean_along_z(r, z, q, conc, species) -> dict:
    """Collapse dofs that share z into one radial mean. Weight r (axisymmetric)."""
    if len(z) == 0:
        return {
            "z_m": [],
            "C_mol_m3": {sp: [] for sp in species},
            "q_mol_kg": [],
        }
    z_sorted = np.asarray(z, dtype=float)
    z_key = np.round(z_sorted, decimals=12)
    order = np.argsort(z_key, kind="mergesort")
    z_key = z_key[order]
    z_sorted = z_sorted[order]
    r = np.asarray(r, dtype=float)[order]
    q = np.asarray(q, dtype=float)[order]
    conc = {sp: np.asarray(conc[sp], dtype=float)[order] for sp in species}
    cuts = np.flatnonzero(np.diff(z_key)) + 1
    spans = np.split(np.arange(len(z_key)), cuts)
    z_out: list[float] = []
    q_out: list[float] = []
    c_out = {sp: [] for sp in species}
    for idx in spans:
        w = np.maximum(r[idx], 1e-12)
        w = w / np.sum(w)
        z_out.append(float(np.dot(z_sorted[idx], w)))
        q_out.append(float(np.dot(q[idx], w)))
        for sp in species:
            c_out[sp].append(float(np.dot(conc[sp][idx], w)))
    return {"z_m": z_out, "C_mol_m3": c_out, "q_mol_kg": q_out}
