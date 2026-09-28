"""
ADR framework: separate problem formulation (Model) from numerical solution.

MVC mapping:
  Model      — AbstractProblem + Advection/Diffusion/Reaction + BC/IC
  View       — plots, JSON export, balance reports
  Controller — SimulationController (CLI / web entry)

Abstract Factory:
  ProblemFactory subclasses build a coherent product family per physical system.
"""

from __future__ import annotations

from adr.model.problem import AbstractProblem, SpeciesSpace
from adr.model.state import FieldStateHistory
from adr.factories.registry import get_factory, list_factories

__all__ = [
    "AbstractProblem",
    "SpeciesSpace",
    "FieldStateHistory",
    "get_factory",
    "list_factories",
]
