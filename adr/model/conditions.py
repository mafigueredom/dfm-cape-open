"""Boundary and initial condition interfaces."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class InletData:
    """Inlet concentrations and optional volumetric flow at time t."""

    C_mol_m3: dict[str, float]
    Q_m3_s: float | None = None
    phase_id: str | None = None
    adsorption_active: bool = True
    methanation_active: bool = True


@dataclass(frozen=True)
class OutletData:
    """Outlet mass BC description (natural Neumann = zero diffusive flux)."""

    kind: str = "zero_diffusive_flux"


@dataclass(frozen=True)
class WallData:
    """Wall mass BC description."""

    kind: str = "zero_diffusive_flux"


class BoundaryConditionSet(ABC):
    """Time-dependent / facet BCs for the ADR problem."""

    @abstractmethod
    def inlet(self, t: float, dt: float = 0.0) -> InletData:
        ...

    @abstractmethod
    def outlet(self) -> OutletData:
        ...

    @abstractmethod
    def wall(self) -> WallData:
        ...

    def inlet_bc_mode(self) -> str:
        """``danckwerts`` | ``dirichlet``."""
        return "danckwerts"

    def apply_to_context(self, ctx: Any, t: float, dt: float) -> InletData:
        """
        Push inlet schedule into the numerical context (velocity, D_z, C_in).

        Default: return inlet data only; DFM overrides to drive FEM fields.
        """
        return self.inlet(t, dt)


class InitialConditionSet(ABC):
    """Initial gas / solid state."""

    @abstractmethod
    def gas_concentrations(self) -> dict[str, float]:
        ...

    @abstractmethod
    def solid_loading(self) -> float:
        """Initial adsorbed loading q [mol/kg] (0 if no solid phase)."""
