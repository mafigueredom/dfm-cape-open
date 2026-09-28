"""
Methanation / DFM concrete Abstract Factory.

Wraps existing POW_V07 physics (prescribed U0, Wilke/DL, CSC+LDF+rxn, Danckwerts).
"""

from __future__ import annotations

from typing import Any

from adr.factories.base import ProblemFactory
from adr.model.conditions import (
    BoundaryConditionSet,
    InitialConditionSet,
    InletData,
    OutletData,
    WallData,
)
from adr.model.problem import AbstractProblem, AuxPhysicsBundle, SpeciesSpace
from adr.model.terms import (
    AdvectionTerm,
    DiffusionTerm,
    OperatorSplitPolicy,
    ReactionTerm,
)
from dfm_config import DfmConfig
from fenicsx.cycle_scheduler import CycleScheduler
from fenicsx.inlet_sopdt import MultiSpeciesInletSOPDT
from fenicsx.properties import (
    diluted_species_initial_C,
    inlet_bc_mode,
    species_list,
)


class DfmAdvectionTerm(AdvectionTerm):
    """Prescribed interstitial plug flow (or Darcy velocity via runtime)."""

    def update_velocity(self, ctx: Any, t: float) -> None:
        # Velocity is refreshed inside BoundaryConditionSet.apply_to_context
        # together with the inlet schedule (matches legacy reactor coupling).
        return None


class DfmDiffusionTerm(DiffusionTerm):
    """Axial dispersion D_z = α D_m + β d_p U_sup."""

    def update_diffusivity(self, ctx: Any, t: float, *, force: bool = False) -> bool:
        return bool(ctx.update_D_z_fields(force=force))


class DfmReactionTerm(ReactionTerm):
    """CSC + LDF + methanation; Lie-split after transport (case A)."""

    def split_policy(self) -> OperatorSplitPolicy:
        return OperatorSplitPolicy.LIE_TRANSPORT_THEN_REACTION

    def transport_source(self, ctx: Any, species: str) -> Any:
        return ctx.zero_source

    def pointwise_exchange(self, ctx: Any, dt: float) -> None:
        ctx.solid_gas_exchange(dt)


class DfmBoundaryConditions(BoundaryConditionSet):
    def __init__(self, cfg: DfmConfig, species: list[str]):
        self.cfg = cfg
        self.species = species
        self.scheduler = CycleScheduler(cfg)
        self._sopdt = MultiSpeciesInletSOPDT(cfg, self.scheduler, species)
        self._sopdt.reset({sp: 0.0 for sp in species})
        self._last_phase_id: str | None = None
        vel = cfg.gas_phase.get("velocity") or {}
        self._Q_default = float(vel.get("Q_m3_s_at_T_ref") or 0.0)
        if self._Q_default <= 0:
            raise ValueError("gas_phase.velocity.Q_m3_s_at_T_ref required.")

    def inlet_bc_mode(self) -> str:
        return inlet_bc_mode(self.cfg)

    def inlet(self, t: float, dt: float = 0.0) -> InletData:
        phase = self.scheduler.phase_at(t)
        if phase.id != self._last_phase_id:
            self._sopdt.begin_phase(t)
            self._last_phase_id = phase.id
        if dt > 0:
            cin = self._sopdt.step(t, dt)
        else:
            cin = self._sopdt.step(t, 1e-6)
        Q = self.scheduler.Q_m3_s_at(t, self._Q_default)
        return InletData(
            C_mol_m3=dict(cin),
            Q_m3_s=Q,
            phase_id=phase.id,
            adsorption_active=bool(phase.adsorption_active),
            methanation_active=bool(phase.methanation_active),
        )

    def outlet(self) -> OutletData:
        return OutletData()

    def wall(self) -> WallData:
        return WallData()

    def apply_to_context(self, ctx: Any, t: float, dt: float) -> InletData:
        return ctx.apply_inlet_schedule(t, dt, self)


class DfmInitialConditions(InitialConditionSet):
    def __init__(self, cfg: DfmConfig, species: list[str]):
        self.cfg = cfg
        self.species = species

    def gas_concentrations(self) -> dict[str, float]:
        return diluted_species_initial_C(self.cfg, self.species)

    def solid_loading(self) -> float:
        q0 = (self.cfg.solid_kinetics.get("q_CO2") or {}).get("initial")
        return float(q0) if q0 is not None else 0.0


class MethanationDfmFactory(ProblemFactory):
    factory_id = "methanation_dfm"
    label = "Isothermal DFM methanation (POW_V07)"

    def __init__(self, config: DfmConfig):
        super().__init__(config)
        self.cfg: DfmConfig = config
        self._species = species_list(config)

    def create_species_space(self) -> SpeciesSpace:
        return SpeciesSpace(
            names=list(self._species),
            stoich_notes="CO2 + 4 H2 -> CH4 + 2 H2O; CSC+LDF solid",
        )

    def create_advection(self) -> AdvectionTerm:
        return DfmAdvectionTerm()

    def create_diffusion(self) -> DiffusionTerm:
        return DfmDiffusionTerm()

    def create_reaction(self) -> ReactionTerm:
        return DfmReactionTerm()

    def create_boundary_conditions(self) -> BoundaryConditionSet:
        return DfmBoundaryConditions(self.cfg, self._species)

    def create_initial_conditions(self) -> InitialConditionSet:
        return DfmInitialConditions(self.cfg, self._species)

    def create_aux_physics(self) -> AuxPhysicsBundle:
        bcs = self.create_boundary_conditions()
        assert isinstance(bcs, DfmBoundaryConditions)
        return AuxPhysicsBundle(
            scheduler=bcs.scheduler,
            property_provider=None,
            config=self.cfg,
            extras={"Q_default": bcs._Q_default},
        )

    def create_problem(self) -> AbstractProblem:
        # Single BC instance shared with aux.scheduler
        bcs = self.create_boundary_conditions()
        assert isinstance(bcs, DfmBoundaryConditions)
        return AbstractProblem(
            species=self.create_species_space(),
            advection=self.create_advection(),
            diffusion=self.create_diffusion(),
            reaction=self.create_reaction(),
            bcs=bcs,
            ics=self.create_initial_conditions(),
            aux=AuxPhysicsBundle(
                scheduler=bcs.scheduler,
                config=self.cfg,
                extras={"Q_default": bcs._Q_default},
            ),
            factory_id=self.factory_id,
            label=self.label,
        )

    def create_runtime(self, problem: AbstractProblem | None = None) -> Any:
        from adr.numerical.dfm_runtime import DfmFenicsRuntime

        prob = problem if problem is not None else self.create_problem()
        return DfmFenicsRuntime(self.cfg, prob)
