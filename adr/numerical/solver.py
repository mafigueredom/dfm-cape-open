"""
AdrSolver / TimeIntegrator — numerical solution of an AbstractProblem.

For methanation_dfm the FEM substrate is ``DfmFenicsRuntime``; the adaptive
time loop and Lie split live here (extracted from IsothermalAxisymmetricReactor).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from adr.model.problem import AbstractProblem
from adr.model.state import FieldStateHistory
from adr.model.terms import OperatorSplitPolicy
from adr.numerical.assembler import WeakFormAssembler
from fenicsx.balance_diagnostics import (
    accumulate_flux_step,
    compute_balance_summary,
    sample_holdups,
)

SampleCallback = Callable[[float, Any, FieldStateHistory], None]


class TimeIntegrator:
    """One physical time step: BC → ADR → reaction (per split policy) → aux."""

    def __init__(self, problem: AbstractProblem, runtime: Any):
        self.problem = problem
        self.runtime = runtime
        self.assembler = WeakFormAssembler(problem)

    def step(self, t: float, dt: float) -> None:
        rt = self.runtime
        self.assembler.ensure_ready(rt)
        self.problem.bcs.apply_to_context(rt, t, dt)
        self.problem.advection.update_velocity(rt, t)

        policy = self.problem.reaction.split_policy()
        for sp in rt.species:
            src = self.problem.reaction.transport_source(rt, sp)
            rt.species_step(sp, dt, src)

        if policy == OperatorSplitPolicy.LIE_TRANSPORT_THEN_REACTION:
            self.problem.reaction.pointwise_exchange(rt, dt)
        # MONOLITHIC_SOURCE: sources already in ADR step

        rt.update_density()
        rt.continuity_step(dt)


class AdrSolver:
    """
    Adaptive driver over FEM ``TimeIntegrator`` or fermentation ODE runtime.

    ``solve(problem, runtime, t_end, dt)`` returns a ``FieldStateHistory``.
    """

    def solve(
        self,
        problem: AbstractProblem,
        runtime: Any,
        *,
        t_end: float,
        dt: float,
        on_sample: SampleCallback | None = None,
    ) -> FieldStateHistory:
        backend = getattr(runtime, "backend", None)
        if backend == "ode":
            return self._solve_ode(
                problem, runtime, t_end=t_end, dt=dt, on_sample=on_sample
            )
        if backend == "fermentation_fem":
            return self._solve_fermentation_fem(
                problem, runtime, t_end=t_end, dt=dt, on_sample=on_sample
            )
        return self._solve_fem(
            problem, runtime, t_end=t_end, dt=dt, on_sample=on_sample
        )

    def _solve_ode(
        self,
        problem: AbstractProblem,
        runtime: Any,
        *,
        t_end: float,
        dt: float,
        on_sample: SampleCallback | None = None,
    ) -> FieldStateHistory:
        """
        Fermentation Phase A: prefer stiff BDF (scipy); RK4 fallback.

        Stabilization: C2 floor, soft O2 inhibition, RHS/state caps
        (see ``kinetic_orders`` in config).
        """
        cfg = problem.aux.config
        sim = (getattr(cfg, "raw", None) or {}).get("simulation") or {}
        dt_step = float(dt)
        record_every = max(1, int(sim.get("record_every_n_steps", 1)))
        method = str(sim.get("ode_method", "BDF")).strip().upper()
        rtol = float(sim.get("ode_rtol", 1e-6))
        atol = float(sim.get("ode_atol", 1e-9))

        state = FieldStateHistory(
            species=list(runtime.species),
            factory_id=problem.factory_id,
        )
        state.outlet_means = {sp: [] for sp in runtime.species}
        state.inlet_effective = {sp: [] for sp in runtime.species}
        state.F_in_mol = {sp: 0.0 for sp in runtime.species}
        state.F_out_mol = {sp: 0.0 for sp in runtime.species}

        def _record_from_y(t_now: float, y: np.ndarray) -> None:
            state.times.append(float(t_now))
            state.phase_id.append("batch")
            state.Q_m3_s.append(0.0)
            for i, sp in enumerate(runtime.species):
                state.outlet_means[sp].append(float(y[i]))
                state.inlet_effective[sp].append(0.0)

        used = method
        t_final = float(t_end)
        if method == "BDF" and hasattr(runtime, "integrate_bdf"):
            try:
                n_rec = max(2, int(np.ceil(t_end / max(dt_step * record_every, 1e-12))) + 1)
                n_rec = min(n_rec, int(sim.get("ode_max_samples", 5001)))
                t_eval = np.linspace(0.0, float(t_end), n_rec)
                problem.bcs.apply_to_context(runtime, 0.0, 0.0)
                t_out, y_out = runtime.integrate_bdf(
                    float(t_end), t_eval=t_eval, rtol=rtol, atol=atol
                )
                for k, tk in enumerate(t_out):
                    _record_from_y(tk, y_out[:, k])
                    if on_sample is not None:
                        on_sample(float(tk), runtime, state)
                t_final = float(t_out[-1])
                used = "BDF"
            except Exception as exc:  # noqa: BLE001 — fall back to RK4
                print(f"  BDF failed ({exc}); falling back to RK4")
                used = "rk4"
                # reset and RK4 below
                runtime.y = np.array(runtime.params["T0"], dtype=float)
                runtime.t = 0.0
                state.times.clear()
                state.phase_id.clear()
                state.Q_m3_s.clear()
                for sp in runtime.species:
                    state.outlet_means[sp].clear()
                    state.inlet_effective[sp].clear()

        if used == "rk4" or method == "RK4":
            t = 0.0
            _record_from_y(t, runtime.y)
            if on_sample is not None:
                on_sample(t, runtime, state)
            step_i = 0
            while t < t_end - 1e-15:
                dt_use = min(dt_step, t_end - t)
                problem.bcs.apply_to_context(runtime, t, dt_use)
                runtime.step_ode(t, dt_use)
                t += dt_use
                step_i += 1
                if step_i % record_every == 0 or t >= t_end - 1e-15:
                    _record_from_y(t, runtime.y)
                    if on_sample is not None:
                        on_sample(t, runtime, state)
            t_final = t
            used = "rk4"

        state.balance_summary = {
            "schema": "fermentation.phase_a",
            "X2": float(runtime.X2),
            "C2_final": float(runtime.C2()),
            "T18_final": float(runtime.T18_at(t_final)),
            "I_final": float(runtime.light.intensity(t_final)),
            "ode_method": used,
            "note": (
                "stabilized S-system: soft O2 inhibition, C2_floor, "
                f"state/RHS caps; integrator={used}"
            ),
        }
        return state

    def _solve_fermentation_fem(
        self,
        problem: AbstractProblem,
        runtime: Any,
        *,
        t_end: float,
        dt: float,
        on_sample: SampleCallback | None = None,
    ) -> FieldStateHistory:
        """Phase B: fixed-dt Lie ADR (diffusion → S-system), domain-mean history."""
        cfg = problem.aux.config
        sim = (getattr(cfg, "raw", None) or {}).get("simulation") or {}
        record_every = max(1, int(sim.get("record_every_n_steps", 1)))
        integrator = TimeIntegrator(problem, runtime)

        state = FieldStateHistory(
            species=list(runtime.species),
            factory_id=problem.factory_id,
        )
        state.outlet_means = {sp: [] for sp in runtime.species}
        state.inlet_effective = {sp: [] for sp in runtime.species}
        state.F_in_mol = {sp: 0.0 for sp in runtime.species}
        state.F_out_mol = {sp: 0.0 for sp in runtime.species}

        def _record(t_now: float) -> None:
            state.times.append(float(t_now))
            state.phase_id.append(str(getattr(runtime, "_current_phase_id", "batch_adr")))
            state.Q_m3_s.append(0.0)
            for sp in runtime.species:
                state.outlet_means[sp].append(runtime.outlet_mean(runtime.C_n[sp]))
                state.inlet_effective[sp].append(0.0)

        t = 0.0
        _record(t)
        if on_sample is not None:
            on_sample(t, runtime, state)
        step_i = 0
        while t < t_end - 1e-15:
            dt_use = min(float(dt), float(t_end) - t)
            integrator.step(t, dt_use)
            t += dt_use
            step_i += 1
            if step_i % record_every == 0 or t >= t_end - 1e-15:
                _record(t)
                if on_sample is not None:
                    on_sample(t, runtime, state)

        o2 = state.outlet_means.get("T2_O2", [0.0])[-1]
        solar = getattr(runtime.light, "solar", None)
        solar_diag: dict[str, float] = {}
        if solar is not None:
            solar_diag = {
                "I_oh_final": float(solar.I_oh(t)),
                "I0_final": float(solar.I0(t)),
                "I_oh_max_48h": float(
                    max(solar.I_oh(ti) for ti in np.linspace(0.0, 48 * 3600.0, 97))
                ),
            }
        state.balance_summary = {
            "schema": "fermentation.phase_b",
            "X2": float(runtime.X2),
            "R_m": float(runtime.R_m),
            "L_m": float(runtime.L_m),
            "O2_mean_final": float(o2),
            "T18_wall": float(runtime.light.T18_path(t, 0.0)),
            "T18_center": float(runtime.light.T18_path(t, runtime.R_m)),
            "I_final": float(runtime.light.intensity(t)),
            "n_dofs": int(runtime.C_n[runtime.species[0]].x.array.shape[0]),
            **solar_diag,
            "note": (
                "spatial ADR U=0: diffusion + Lie S-system; "
                "I0(t) from Liu–Jordan I_oh; light_path=radial_wall (P_dis=R-r)"
            ),
        }
        return state

    def _solve_fem(
        self,
        problem: AbstractProblem,
        runtime: Any,
        *,
        t_end: float,
        dt: float,
        on_sample: SampleCallback | None = None,
    ) -> FieldStateHistory:
        cfg = problem.aux.config
        if cfg is None:
            raise ValueError("AbstractProblem.aux.config required for AdrSolver")

        integrator = TimeIntegrator(problem, runtime)
        sim = cfg.raw.get("simulation") or {}
        dt_min = float(dt)
        dt_max_req = float(sim.get("dt_max_s") or dt_min)
        dt_max = max(min(dt_max_req, 0.1), dt_min)
        adaptive = dt_max > dt_min * (1.0 + 1e-9)
        tol_lo = float(sim.get("adapt_tol_lo", 0.001))
        tol_hi = float(sim.get("adapt_tol_hi", 0.005))
        c_floor = max(
            0.01
            * max(
                (
                    float(np.max(np.abs(runtime.C_n[sp].x.array)))
                    for sp in runtime.species
                ),
                default=1.0,
            ),
            1e-9,
        )

        scheduler = problem.aux.scheduler
        event_times: list[float] = []
        if scheduler is not None and hasattr(scheduler, "boundary_times"):
            event_times = list(scheduler.boundary_times(t_end=t_end))
        sopdt = (cfg.gas_phase.get("boundaries_mass") or {}).get("z_equals_0") or {}
        sopdt = sopdt.get("sopdt") or {}
        theta_p = float(sopdt.get("theta_p", 0.0)) if sopdt.get("enabled") else 0.0
        if theta_p > 0:
            event_times.append(theta_p)
        event_times = sorted(set(event_times))
        event_pad = float(sim.get("adapt_event_pad_s", 30.0))
        boundaries = sorted(set(event_times) | {float(t_end)})

        Q_default = float(problem.aux.extras.get("Q_default") or 0.0)
        state = FieldStateHistory(
            species=list(runtime.species),
            factory_id=problem.factory_id,
        )
        state.outlet_means = {sp: [] for sp in runtime.species}
        state.inlet_effective = {sp: [] for sp in runtime.species}
        state.F_in_mol = {sp: 0.0 for sp in runtime.species}
        state.F_out_mol = {sp: 0.0 for sp in runtime.species}
        holdup_keys = list(runtime.species) + ["CO2_ads"]
        if runtime.balance_on:
            state.holdup_mol = {k: [] for k in holdup_keys}

        t = 0.0
        Q_now = scheduler.Q_m3_s_at(t, Q_default)
        state.times.append(t)
        state.phase_id.append(runtime._current_phase_id)
        state.Q_m3_s.append(Q_now)
        for sp in runtime.species:
            state.outlet_means[sp].append(runtime.outlet_mean(runtime.C_n[sp]))
            state.inlet_effective[sp].append(runtime.C_in[sp])
        holdup_0: dict[str, float] | None = None
        if runtime.balance_on and runtime.holdup is not None:
            h0 = sample_holdups(
                runtime.holdup, runtime.C_n, runtime.q_n, runtime.species
            )
            holdup_0 = dict(h0)
            for k in holdup_keys:
                state.holdup_mol[k].append(float(h0[k]))

        if on_sample is not None:
            on_sample(t, runtime, state)

        dt_cur = dt_min
        prev = (
            {sp: runtime.C_n[sp].x.array.copy() for sp in runtime.species}
            if adaptive
            else {}
        )
        cin_prev = {sp: float(runtime.C_in[sp]) for sp in runtime.species}
        Q_prev = Q_now
        Cin_prev = {sp: float(runtime.C_in[sp]) for sp in runtime.species}
        Cout_prev = {sp: float(state.outlet_means[sp][-1]) for sp in runtime.species}
        step_i = 0

        while t < t_end - 1e-9:
            t_next = next((b for b in boundaries if b > t + 1e-9), t_end)
            near_event = any(
                abs(t - te) < event_pad or 0 <= (t - te) < event_pad
                for te in event_times
            )
            if near_event:
                dt_cur = dt_min
            dt_step = min(dt_cur, t_next - t)
            if adaptive:
                for sp in runtime.species:
                    prev[sp][:] = runtime.C_n[sp].x.array

            integrator.step(t, dt_step)
            t += dt_step
            step_i += 1

            Q_now = scheduler.Q_m3_s_at(t, Q_default)
            Cout_now = {
                sp: runtime.outlet_mean(runtime.C_n[sp]) for sp in runtime.species
            }
            Cin_now = {sp: float(runtime.C_in[sp]) for sp in runtime.species}
            if runtime.balance_on:
                accumulate_flux_step(
                    state.F_in_mol,
                    species=runtime.species,
                    Q0=Q_prev,
                    Q1=Q_now,
                    C0=Cin_prev,
                    C1=Cin_now,
                    dt=dt_step,
                )
                accumulate_flux_step(
                    state.F_out_mol,
                    species=runtime.species,
                    Q0=Q_prev,
                    Q1=Q_now,
                    C0=Cout_prev,
                    C1=Cout_now,
                    dt=dt_step,
                )

            state.times.append(t)
            state.phase_id.append(scheduler.phase_at(t).id)
            state.Q_m3_s.append(Q_now)
            for sp in runtime.species:
                state.outlet_means[sp].append(Cout_now[sp])
                state.inlet_effective[sp].append(Cin_now[sp])

            if on_sample is not None:
                on_sample(t, runtime, state)

            sample_holdup = runtime.balance_on and (
                step_i % runtime.balance_every == 0 or t >= t_end - 1e-9
            )
            if sample_holdup and runtime.holdup is not None:
                h = sample_holdups(
                    runtime.holdup, runtime.C_n, runtime.q_n, runtime.species
                )
                for k in holdup_keys:
                    state.holdup_mol[k].append(float(h[k]))
            elif runtime.balance_on:
                for k in holdup_keys:
                    state.holdup_mol[k].append(state.holdup_mol[k][-1])

            Q_prev = Q_now
            Cin_prev = Cin_now
            Cout_prev = Cout_now

            if adaptive:
                if abs(t - t_next) < 1e-9 and t < t_end - 1e-9:
                    dt_cur = dt_min
                    cin_prev = {sp: float(runtime.C_in[sp]) for sp in runtime.species}
                    continue
                metric = 0.0
                for sp in runtime.species:
                    arr = runtime.C_n[sp].x.array
                    scale = max(float(np.max(np.abs(arr))), c_floor)
                    d = float(np.max(np.abs(arr - prev[sp])))
                    metric = max(metric, d / scale)
                cin_chg = max(
                    abs(float(runtime.C_in[sp]) - cin_prev[sp])
                    / max(abs(float(runtime.C_in[sp])), c_floor)
                    for sp in runtime.species
                )
                cin_prev = {sp: float(runtime.C_in[sp]) for sp in runtime.species}
                if near_event or cin_chg > tol_hi:
                    dt_cur = dt_min
                elif metric > tol_hi:
                    dt_cur = max(dt_min, 0.5 * dt_cur)
                elif metric < tol_lo:
                    dt_cur = min(dt_max, 2.0 * dt_cur)

        if runtime.balance_on and holdup_0 is not None and state.holdup_mol:
            holdup_f = {k: float(state.holdup_mol[k][-1]) for k in holdup_keys}
            state.balance_summary = compute_balance_summary(
                species=runtime.species,
                holdup_0=holdup_0,
                holdup_f=holdup_f,
                F_in=state.F_in_mol,
                F_out=state.F_out_mol,
            )
        return state
