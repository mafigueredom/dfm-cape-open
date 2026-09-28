"""Piecewise-constant inlet / phase flags from simulation.cycle in dfm_config."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from dfm_config import DfmConfig
from fenicsx.properties import species_list


def n_cycles_from_sim(sim: dict[str, Any] | None) -> int:
    sim = sim or {}
    try:
        n = int(sim.get("n_cycles") or 1)
    except (TypeError, ValueError):
        n = 1
    return max(1, n)


def cycle_period_s_from_sim(sim: dict[str, Any] | None) -> float:
    """Period of one cycle from phase windows (max t_end − min t_start)."""
    sim = sim or {}
    phases = (sim.get("cycle") or {}).get("phases") or []
    usable: list[tuple[float, float]] = []
    for p in phases:
        t0 = float(p.get("t_start_s", 0.0))
        t1 = float(p.get("t_end_s", 0.0))
        if t1 > t0 + 1e-12:
            usable.append((t0, t1))
    if not usable:
        return 0.0
    return max(t1 for _, t1 in usable) - min(t0 for t0, _ in usable)


def resolved_t_end_s(sim: dict[str, Any] | None) -> float:
    """FEM horizon: n_cycles × period when n_cycles > 1, else t_end_s."""
    sim = sim or {}
    t_end = float(sim.get("t_end_s") or 0.0)
    n = n_cycles_from_sim(sim)
    period = cycle_period_s_from_sim(sim)
    if n > 1 and period > 0.0:
        return float(n) * period
    return t_end


@dataclass(frozen=True)
class Phase:
    id: str
    t_start_s: float
    t_end_s: float
    inlet_C_mol_m3: dict[str, float]
    methanation_active: bool = True
    adsorption_active: bool = True
    Q_m3_s: float | None = None


class CycleScheduler:
    def __init__(self, cfg: DfmConfig):
        self.cfg = cfg
        self.species = species_list(cfg)
        sim = cfg.raw.get("simulation") or {}
        self.n_cycles = n_cycles_from_sim(sim)
        cycle = sim.get("cycle")
        if not cycle or not cycle.get("phases"):
            fallback = sim.get("inlet_C_mol_m3") or {}
            cin = {sp: float(fallback[sp]) for sp in self.species}
            self.phases = [
                Phase(
                    id="steady_inlet",
                    t_start_s=0.0,
                    t_end_s=float(sim.get("t_end_s", 1e9)),
                    inlet_C_mol_m3=cin,
                    methanation_active=True,
                    adsorption_active=True,
                )
            ]
            self._static = True
            self._t0 = 0.0
            self._period = float(self.phases[0].t_end_s)
            return
        self._static = False
        self.phases = []
        for raw in cycle["phases"]:
            t0 = float(raw["t_start_s"])
            t1 = float(raw["t_end_s"])
            if t1 <= t0 + 1e-12:
                continue
            cin_raw = raw.get("inlet_C_mol_m3") or {}
            cin = {sp: float(cin_raw[sp]) for sp in self.species}
            self.phases.append(
                Phase(
                    id=str(raw.get("id", "phase")),
                    t_start_s=t0,
                    t_end_s=t1,
                    inlet_C_mol_m3=cin,
                    methanation_active=bool(raw.get("methanation_active", True)),
                    adsorption_active=bool(raw.get("adsorption_active", True)),
                    Q_m3_s=(
                        float(raw["Q_m3_s"])
                        if raw.get("Q_m3_s") is not None
                        else None
                    ),
                )
            )
        if not self.phases:
            raise ValueError("simulation.cycle.phases has no positive-duration steps")
        self.phases.sort(key=lambda p: p.t_start_s)
        self._t0 = float(self.phases[0].t_start_s)
        self._period = float(self.phases[-1].t_end_s) - self._t0
        self._validate()

    def _validate(self) -> None:
        for p in self.phases:
            if p.t_end_s <= p.t_start_s:
                raise ValueError(f"cycle phase {p.id}: t_end must exceed t_start")
        for i in range(len(self.phases) - 1):
            if self.phases[i].t_end_s > self.phases[i + 1].t_start_s + 1e-9:
                raise ValueError("cycle phases overlap in time")

    def t_cycle_s(self) -> float:
        return float(self._period)

    def _tau(self, t: float) -> float:
        """Map absolute time into the first-cycle window when repeating."""
        if self.n_cycles <= 1 or self._period <= 0.0:
            return t
        rel = float(t) - self._t0
        if rel < 0.0:
            return t
        tau = rel % self._period
        return self._t0 + tau

    def phase_at(self, t: float) -> Phase:
        t_loc = self._tau(t)
        for p in self.phases:
            if p.t_start_s <= t_loc < p.t_end_s - 1e-12:
                return p
        return self.phases[-1]

    def inlet_at(self, t: float) -> dict[str, float]:
        if t < 0.0:
            return {sp: 0.0 for sp in self.species}
        return dict(self.phase_at(t).inlet_C_mol_m3)

    def methanation_active_at(self, t: float) -> bool:
        return bool(self.phase_at(t).methanation_active)

    def adsorption_active_at(self, t: float) -> bool:
        return bool(self.phase_at(t).adsorption_active)

    def Q_m3_s_at(self, t: float, default_Q: float) -> float:
        q = self.phase_at(t).Q_m3_s
        return float(q) if q is not None and q > 0 else default_Q

    def boundary_times(self, t_end: float | None = None) -> list[float]:
        """Phase-switch times in [0, t_end) for adaptive stepping / plots."""
        T = self.t_cycle_s()
        n = self.n_cycles if T > 0.0 else 1
        out: list[float] = []
        for k in range(max(1, n)):
            off = k * T if T > 0.0 else 0.0
            if k > 0:
                out.append(self._t0 + off)
            for p in self.phases:
                ts = p.t_start_s + off
                if ts > 0.0:
                    out.append(ts)
        uniq = sorted(set(out))
        if t_end is None:
            return uniq
        te = float(t_end)
        return [x for x in uniq if x < te - 1e-12]

    def phase_labels(self) -> list[tuple[float, str]]:
        T = self.t_cycle_s()
        n = self.n_cycles if T > 0.0 else 1
        out: list[tuple[float, str]] = []
        for k in range(max(1, n)):
            off = k * T if T > 0.0 else 0.0
            suffix = "" if n <= 1 else f"#{k + 1}"
            for p in self.phases:
                out.append((p.t_start_s + off, f"{p.id}{suffix}"))
        return out
