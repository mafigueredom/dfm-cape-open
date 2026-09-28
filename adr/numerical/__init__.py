"""Numerical solution layer (dolfinx assemblers, time integrators)."""

from __future__ import annotations

from adr.model.terms import OperatorSplitPolicy
from adr.numerical.solver import AdrSolver

__all__ = ["AdrSolver", "OperatorSplitPolicy"]
