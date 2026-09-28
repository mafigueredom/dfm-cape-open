"""Time-history state produced by the numerical solver (View-consumable)."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class FieldStateHistory:
    """
    Observable time series from an ADR solve.

    Compatible with the former ``ReactorState`` used by export/plots.
    """

    times: list[float] = field(default_factory=list)
    outlet_means: dict[str, list[float]] = field(default_factory=dict)
    inlet_effective: dict[str, list[float]] = field(default_factory=dict)
    phase_id: list[str] = field(default_factory=list)
    Q_m3_s: list[float] = field(default_factory=list)
    holdup_mol: dict[str, list[float]] = field(default_factory=dict)
    F_in_mol: dict[str, float] = field(default_factory=dict)
    F_out_mol: dict[str, float] = field(default_factory=dict)
    balance_summary: dict | None = None
    species: list[str] = field(default_factory=list)
    factory_id: str = ""


# Backward-compatible alias used across fenicsx export / plots
ReactorState = FieldStateHistory
