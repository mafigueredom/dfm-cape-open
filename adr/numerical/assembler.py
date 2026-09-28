"""
Weak-form assembly helpers for ADR transport.

DFM currently builds cached ``SpeciesTransportSolver`` instances inside
``DfmFenicsRuntime`` (same forms as ``fenicsx/diluted_species.py``). This
module documents the assembler role for future non-DFM backends.
"""

from __future__ import annotations

from typing import Any

from adr.model.problem import AbstractProblem
from adr.model.terms import OperatorSplitPolicy


class WeakFormAssembler:
    """Bind an AbstractProblem to backend transport solvers."""

    def __init__(self, problem: AbstractProblem):
        self.problem = problem

    def transport_uses_zero_kinetic_source(self) -> bool:
        return (
            self.problem.reaction.split_policy()
            == OperatorSplitPolicy.LIE_TRANSPORT_THEN_REACTION
        )

    def ensure_ready(self, runtime: Any) -> None:
        """No-op for DFM (solvers constructed in runtime.__init__)."""
        if not hasattr(runtime, "transport_solvers"):
            raise RuntimeError("runtime missing transport_solvers")
