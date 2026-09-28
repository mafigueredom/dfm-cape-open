"""Abstract ADR problem: assembled formulation products from a ProblemFactory."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from adr.model.conditions import BoundaryConditionSet, InitialConditionSet
from adr.model.terms import AdvectionTerm, DiffusionTerm, ReactionTerm


@dataclass(frozen=True)
class SpeciesSpace:
    """Gas (and optional solid) species participating in the ADR system."""

    names: list[str]
    stoich_notes: str = ""


@dataclass
class AuxPhysicsBundle:
    """
    Extra physics attached to a problem family (scheduler, solid ODE hooks, cfg).

    Kept intentionally loose so SystemB/C can store different auxiliaries.
    """

    scheduler: Any = None
    property_provider: Any = None
    config: Any = None
    extras: dict[str, Any] = field(default_factory=dict)


@dataclass
class AbstractProblem:
    """Complete ADR formulation for one physical system."""

    species: SpeciesSpace
    advection: AdvectionTerm
    diffusion: DiffusionTerm
    reaction: ReactionTerm
    bcs: BoundaryConditionSet
    ics: InitialConditionSet
    aux: AuxPhysicsBundle
    factory_id: str = ""
    label: str = ""
