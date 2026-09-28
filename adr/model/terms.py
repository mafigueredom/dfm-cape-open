"""ADR term interfaces: advection, diffusion, reaction."""

from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from typing import Any


class OperatorSplitPolicy(str, Enum):
    """How the numerical layer applies ReactionTerm relative to transport."""

    MONOLITHIC_SOURCE = "monolithic_source"
    LIE_TRANSPORT_THEN_REACTION = "lie_transport_then_reaction"


class AdvectionTerm(ABC):
    """Formulation of the advective flux / velocity field."""

    @abstractmethod
    def update_velocity(self, ctx: Any, t: float) -> None:
        """Refresh velocity on ``ctx`` for time ``t`` (may be a no-op)."""

    def contribute_weak_form(self, trial: Any, test: Any, ctx: Any) -> Any | None:
        """Optional UFL contribution; default forms may already include advection."""
        return None


class DiffusionTerm(ABC):
    """Formulation of dispersive / diffusive flux."""

    @abstractmethod
    def update_diffusivity(self, ctx: Any, t: float, *, force: bool = False) -> bool:
        """
        Refresh diffusivity fields on ``ctx``.

        Returns True if coefficients changed enough to require matrix rebuild.
        """

    def contribute_weak_form(self, trial: Any, test: Any, ctx: Any) -> Any | None:
        return None


class ReactionTerm(ABC):
    """
    Formulation of kinetic / exchange sources R_i.

    The term belongs to the Model. ``split_policy`` tells the numerical layer
    whether to put R_i in the ADR weak form or apply a Lie split after transport.
    """

    @abstractmethod
    def split_policy(self) -> OperatorSplitPolicy:
        ...

    @abstractmethod
    def transport_source(self, ctx: Any, species: str) -> Any:
        """Source array/field for species in the ADR step (zeros under Lie split)."""

    @abstractmethod
    def pointwise_exchange(self, ctx: Any, dt: float) -> None:
        """Lie-split solid–gas exchange (no-op if policy is monolithic)."""
