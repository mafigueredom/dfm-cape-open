"""Abstract Factory for ADR problem product families."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from adr.model.conditions import BoundaryConditionSet, InitialConditionSet
from adr.model.problem import AbstractProblem, AuxPhysicsBundle, SpeciesSpace
from adr.model.terms import AdvectionTerm, DiffusionTerm, ReactionTerm


class ProblemFactory(ABC):
    """
    Creates a coherent set of ADR products for one physical system.

    Subclasses must not mix terms across systems (e.g. DFM diffusion + fermentation reaction).
    """

    factory_id: str = ""
    label: str = ""

    def __init__(self, config: Any):
        self.config = config

    @abstractmethod
    def create_advection(self) -> AdvectionTerm:
        ...

    @abstractmethod
    def create_diffusion(self) -> DiffusionTerm:
        ...

    @abstractmethod
    def create_reaction(self) -> ReactionTerm:
        ...

    @abstractmethod
    def create_boundary_conditions(self) -> BoundaryConditionSet:
        ...

    @abstractmethod
    def create_initial_conditions(self) -> InitialConditionSet:
        ...

    @abstractmethod
    def create_species_space(self) -> SpeciesSpace:
        ...

    @abstractmethod
    def create_aux_physics(self) -> AuxPhysicsBundle:
        ...

    def create_problem(self) -> AbstractProblem:
        return AbstractProblem(
            species=self.create_species_space(),
            advection=self.create_advection(),
            diffusion=self.create_diffusion(),
            reaction=self.create_reaction(),
            bcs=self.create_boundary_conditions(),
            ics=self.create_initial_conditions(),
            aux=self.create_aux_physics(),
            factory_id=self.factory_id,
            label=self.label,
        )

    def create_runtime(self, problem: AbstractProblem | None = None) -> Any:
        """
        Optional FEM/runtime backend for this factory family.

        Default: None (solver may reject). DFM overrides.
        """
        return None
