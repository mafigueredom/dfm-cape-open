"""SimulationController — MVC Controller for ADR runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from mpi4py import MPI

from adr.factories.registry import get_factory
from adr.model.state import FieldStateHistory
from adr.numerical.solver import AdrSolver
from adr.view.export import ResultExportView


@dataclass
class SimulationResult:
    state: FieldStateHistory
    species: list[str]
    cfg: Any
    factory_id: str
    config_path: Path | None = None
    phase_markers: list[tuple[float, str]] | None = None
    adsorption_end_s: float = 300.0


def _time_stepping(cfg: Any) -> tuple[float, float]:
    from fenicsx.cycle_scheduler import resolved_t_end_s

    sim = (getattr(cfg, "raw", None) or {}).get("simulation") or {}
    dt = float(sim.get("dt_s", 0.01))
    t_end = resolved_t_end_s(sim)
    if t_end <= 0 or dt <= 0:
        raise ValueError("simulation.t_end_s (or n_cycles × period) and dt_s must be positive")
    return t_end, dt


class SimulationController:
    """
    Load config → Abstract Factory → AdrSolver → optional View.

    Entry points (CLI / web) should call ``run`` rather than constructing
    reactors or solvers directly.
    """

    def __init__(self, *, factory_id: str | None = None):
        self.factory_id = factory_id
        self.solver = AdrSolver()

    def run(
        self,
        config: Any,
        *,
        factory_id: str | None = None,
        no_balance: bool = False,
        t_end: float | None = None,
        dt: float | None = None,
        export_json: Path | str | None = None,
        plot: bool = False,
        results_dir: Path | str | None = None,
        print_summary: bool = True,
        on_sample: Any = None,
    ) -> SimulationResult:
        if isinstance(config, (str, Path)):
            config_path = Path(config)
            from adr.fermentation.config import load_any_config

            cfg = load_any_config(config_path)
        else:
            cfg = config
            config_path = None

        if no_balance:
            sim = cfg.raw.setdefault("simulation", {})
            bal = sim.get("balance")
            if isinstance(bal, dict):
                bal["enabled"] = False
            else:
                sim["balance"] = {"enabled": False}

        fid = factory_id if factory_id is not None else self.factory_id
        factory = get_factory(fid, cfg)
        problem = factory.create_problem()
        runtime = factory.create_runtime(problem)

        t_end_run, dt_run = _time_stepping(cfg)
        if t_end is not None:
            t_end_run = float(t_end)
        if dt is not None:
            dt_run = float(dt)

        state = self.solver.solve(
            problem, runtime, t_end=t_end_run, dt=dt_run, on_sample=on_sample
        )

        scheduler = problem.aux.scheduler
        phase_markers = scheduler.phase_labels() if scheduler is not None else []
        adsorption_end_s = (
            float(scheduler.phases[0].t_end_s)
            if scheduler and getattr(scheduler, "phases", None)
            else 300.0
        )
        sim = cfg.raw.get("simulation") or {}
        comsol_steps = sim.get("comsol_steps") or {}
        if comsol_steps.get("step1_CO2_OFF_location_s"):
            adsorption_end_s = float(comsol_steps["step1_CO2_OFF_location_s"])
        elif (sim.get("cycle") or {}).get("adsorption_end_s"):
            adsorption_end_s = float(sim["cycle"]["adsorption_end_s"])

        result = SimulationResult(
            state=state,
            species=list(runtime.species),
            cfg=cfg,
            factory_id=problem.factory_id,
            config_path=config_path,
            phase_markers=phase_markers,
            adsorption_end_s=adsorption_end_s,
        )

        if MPI.COMM_WORLD.rank == 0 and print_summary:
            self._print_summary(result)

        if MPI.COMM_WORLD.rank == 0:
            view_export = ResultExportView()
            if export_json is not None:
                out = view_export.export_json(
                    state,
                    path=export_json,
                    species=result.species,
                    config_path=config_path,
                    cfg=cfg,
                )
                print(f"Exported run: {out}")
            view_export.print_balance(state)

            if plot:
                if MPI.COMM_WORLD.size != 1:
                    print("Skipping plot: use mpirun -n 1 for --plot")
                else:
                    rdir = (
                        Path(results_dir)
                        if results_dir
                        else Path(__file__).resolve().parents[2] / "results"
                    )
                    if result.factory_id == "fermentation":
                        from adr.view.fermentation_plots import plot_fermentation

                        paths = plot_fermentation(result, rdir, runtime=runtime)
                    else:
                        from adr.view.plots import ResultPlotView

                        paths = ResultPlotView().render(result, results_dir=rdir)
                    print(f"Wrote {len(paths)} figures under {rdir}/")
                    for p in paths:
                        print(f"  {p.name}")

        return result

    @staticmethod
    def _print_summary(result: SimulationResult) -> None:
        state = result.state
        print(
            f"Completed {len(state.times) - 1} recorded samples, "
            f"t_end={state.times[-1]:.4g} s  [factory={result.factory_id}]"
        )
        if result.phase_markers and len(result.phase_markers) > 1:
            print(
                "  Cycle phases:",
                ", ".join(f"{lab} @{t}s" for t, lab in result.phase_markers),
            )
        highlight = {"fermentation"}
        if result.factory_id in highlight:
            keys = [
                "T2_O2",
                "T5_starch",
                "T11_H2",
                "T12_ATP",
                "T16_PSII",
            ]
            for sp in keys:
                if sp not in state.outlet_means:
                    continue
                y0 = state.outlet_means[sp][0]
                yf = state.outlet_means[sp][-1]
                print(f"  {sp}: {y0:.4g} -> {yf:.4g}")
            if state.balance_summary:
                bs = state.balance_summary
                if bs.get("schema") == "fermentation.phase_b":
                    print(
                        f"  ADR: n_dofs={bs.get('n_dofs')}  "
                        f"T18_wall={bs.get('T18_wall', float('nan')):.4g}  "
                        f"T18_center={bs.get('T18_center', float('nan')):.4g}  "
                        f"O2_mean={bs.get('O2_mean_final', float('nan')):.4g}"
                    )
                    if "I_oh_max_48h" in bs:
                        print(
                            f"  solar: I_oh_max(48h)={bs['I_oh_max_48h']:.4g} W/m^2  "
                            f"I0(t_end)={bs.get('I0_final', float('nan')):.4g}"
                        )
                else:
                    print(
                        f"  light: I={bs.get('I_final', float('nan')):.4g}  "
                        f"T18={bs.get('T18_final', float('nan')):.4g}  "
                        f"C2={bs.get('C2_final', float('nan')):.4g}"
                    )
            return
        for sp in result.species:
            y0 = state.outlet_means[sp][0]
            yf = state.outlet_means[sp][-1]
            print(f"  outlet mean {sp}: {y0:.4g} -> {yf:.4g} mol/m^3")
