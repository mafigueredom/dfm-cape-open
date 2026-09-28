"""Problem formulation layer (physics model, backend-agnostic interfaces)."""

from __future__ import annotations

from adr.model.conditions import BoundaryConditionSet, InitialConditionSet, InletData
from adr.model.problem import AbstractProblem, AuxPhysicsBundle, SpeciesSpace
from adr.model.state import FieldStateHistory
from adr.model.terms import (
    AdvectionTerm,
    DiffusionTerm,
    OperatorSplitPolicy,
    ReactionTerm,
)

__all__ = [
    "AbstractProblem",
    "AdvectionTerm",
    "AuxPhysicsBundle",
    "BoundaryConditionSet",
    "DiffusionTerm",
    "FieldStateHistory",
    "InitialConditionSet",
    "InletData",
    "OperatorSplitPolicy",
    "ReactionTerm",
    "SpeciesSpace",
]
