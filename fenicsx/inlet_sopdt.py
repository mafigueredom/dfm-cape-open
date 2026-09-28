"""Inlet concentration dynamics: per-species SOPDT + dead time on phase setpoints."""

from __future__ import annotations

import numpy as np

from dfm_config import DfmConfig

from fenicsx.cycle_scheduler import CycleScheduler


def _sopdt_block(cfg: DfmConfig) -> dict | None:
    bc = (cfg.gas_phase.get("boundaries_mass") or {}).get("z_equals_0") or {}
    sopdt = bc.get("sopdt") or {}
    if not sopdt.get("enabled", False):
        return None
    return sopdt


def _sopdt_params(cfg: DfmConfig) -> dict | None:
    sopdt = _sopdt_block(cfg)
    if sopdt is None:
        return None
    return {
        "Kp": float(sopdt.get("Kp", 1.0)),
        "tau_s": float(sopdt.get("tau_s", 1.0)),
        "zeta": float(sopdt.get("zeta", 1.0)),
        "theta_p": float(sopdt.get("theta_p", 0.0)),
        "per_species": sopdt.get("per_species") or {},
    }


def _row_for_direction(row: dict, rising: bool, fallback: dict) -> dict | None:
    if "rising" in row or "falling" in row:
        key = "rising" if rising else "falling"
        block = row.get(key)
        if not isinstance(block, dict):
            return None
        return block
    if any(k in row for k in ("tau_s", "zeta", "theta_p", "Kp")):
        return row
    return None


def _sopdt_quadruple_for_species(
    sp: str,
    params: dict,
    *,
    rising: bool,
) -> tuple[float, float, float, float] | None:
    """Return (Kp, tau_s, zeta, theta_p) for species and signal direction."""
    row = (params.get("per_species") or {}).get(sp)
    if row is None:
        return None
    block = _row_for_direction(row, rising, params)
    if block is None:
        return None
    return (
        float(block.get("Kp", params["Kp"])),
        float(block.get("tau_s", params["tau_s"])),
        float(block.get("zeta", params["zeta"])),
        float(block.get("theta_p", params["theta_p"])),
    )


class MultiSpeciesInletSOPDT:
    """
    Second-order lag with per-species dead time on scheduler setpoints.

    d²y/dt² + (2ζ/τ) dy/dt + (1/τ²) y = (Kp/τ²) u(t)

    Config may define ``per_species[sp].rising`` and ``per_species[sp].falling``
    (lab SOPDT fits for valve ON/OFF). Species without ``per_species`` pass through
    the raw cycle setpoint (no lag).

    Dead time θ applies on **rising** setpoints. **Falling** setpoints use the
    immediate phase command; decay uses falling (τ, ζ, Kp).
    """

    def __init__(self, cfg: DfmConfig, scheduler: CycleScheduler, species: list[str]):
        self.scheduler = scheduler
        self.species = species
        self.params = _sopdt_params(cfg)
        self.enabled = self.params is not None
        self.y = {sp: 0.0 for sp in species}
        self.dy = {sp: 0.0 for sp in species}
        self._phase_start_t = 0.0
        self._phase_local_theta = False
        self._phase_prev_inlet: dict[str, float] = {sp: 0.0 for sp in species}

    def begin_phase(self, t: float) -> None:
        phase = self.scheduler.phase_at(t)
        self._phase_local_theta = bool(phase.methanation_active)
        if self._phase_local_theta:
            self._phase_start_t = float(t)
            # Hold composition from just before the switch during dead time
            # (inlet_at(t0) is already the NEW phase setpoint).
            t_prev = max(float(t) - 1e-9, 0.0)
            prev = self.scheduler.inlet_at(t_prev)
            self._phase_prev_inlet = {sp: float(prev.get(sp, 0.0)) for sp in self.species}

    def reset(self, y0: dict[str, float] | None = None) -> None:
        if y0:
            for sp in self.species:
                self.y[sp] = float(y0.get(sp, 0.0))
        else:
            for sp in self.species:
                self.y[sp] = 0.0
        for sp in self.species:
            self.dy[sp] = 0.0

    def _delayed_setpoint(self, sp: str, t: float, theta: float) -> float:
        t0 = self._phase_start_t
        if self._phase_local_theta and theta > 0.0 and t - t0 < theta:
            return float(self._phase_prev_inlet.get(sp, 0.0))
        if theta > 0.0:
            return float(self.scheduler.inlet_at(t - theta)[sp])
        return float(self.scheduler.inlet_at(t)[sp])

    def _is_rising(self, sp: str, t: float, theta_rise: float) -> bool:
        immediate = float(self.scheduler.inlet_at(t)[sp])
        delayed = self._delayed_setpoint(sp, t, theta_rise)
        return immediate > delayed + 1e-15

    def _command_at(self, t: float) -> dict[str, float]:
        if not self.enabled:
            return self.scheduler.inlet_at(t)
        immediate = self.scheduler.inlet_at(t)
        out: dict[str, float] = {}
        for sp in self.species:
            quad_rise = _sopdt_quadruple_for_species(sp, self.params, rising=True)
            if quad_rise is None:
                out[sp] = float(immediate[sp])
                continue
            _, _, _, theta_rise = quad_rise
            imm = float(immediate[sp])
            if self._is_rising(sp, t, theta_rise):
                out[sp] = self._delayed_setpoint(sp, t, theta_rise)
            else:
                out[sp] = imm
        return out

    def _euler_step(self, t: float, dt: float) -> None:
        """One forward-Euler step of the SOPDT state at time t with size dt."""
        u = self._command_at(t)
        for sp in self.species:
            quad_rise = _sopdt_quadruple_for_species(sp, self.params, rising=True)
            if quad_rise is None:
                self.y[sp] = float(u[sp])
                self.dy[sp] = 0.0
                continue
            rising = self._is_rising(sp, t, quad_rise[3])
            quad = _sopdt_quadruple_for_species(sp, self.params, rising=rising)
            if quad is None:
                quad = quad_rise
            Kp, tau, zeta, _ = quad
            inv_tau2 = 1.0 / (tau * tau)
            target = Kp * u[sp]
            y = self.y[sp]
            dy = self.dy[sp]
            ddy = inv_tau2 * (target - y) - (2.0 * zeta / tau) * dy
            dy = dy + dt * ddy
            y = y + dt * dy
            y = max(y, 0.0)
            if target > 0.0:
                y = min(y, target)
            self.y[sp] = y
            self.dy[sp] = dy

    def step(self, t: float, dt: float) -> dict[str, float]:
        """
        Advance SOPDT states over [t, t+dt] (caller passes the PDE step start time).

        Substeps with dt_sub ≤ min(0.05, τ/20) so the forward-Euler filter stays
        accurate when the PDE uses a larger adaptive step.
        """
        if not self.enabled:
            return self.scheduler.inlet_at(t)
        if dt <= 0.0:
            return dict(self.y)

        tau_min = None
        for sp in self.species:
            for rising in (True, False):
                quad = _sopdt_quadruple_for_species(sp, self.params, rising=rising)
                if quad is not None:
                    tau_min = float(quad[1]) if tau_min is None else min(tau_min, float(quad[1]))
        dt_cap = 0.05 if tau_min is None else min(0.05, max(tau_min / 20.0, 1e-3))
        n_sub = max(1, int(np.ceil(float(dt) / dt_cap)))
        dt_sub = float(dt) / n_sub
        for k in range(n_sub):
            self._euler_step(t + k * dt_sub, dt_sub)
        return dict(self.y)


def _interp_u(t: np.ndarray, y: np.ndarray, t_query: float, t0: float) -> float:
    if t_query <= t0:
        return 0.0
    return float(np.interp(t_query, t, y))


def apply_sopdt_filter(
    t: np.ndarray,
    y: np.ndarray,
    *,
    Kp: float = 1.0,
    tau_s: float = 1.0,
    zeta: float = 1.0,
    theta_p: float = 0.0,
    dt_sub: float = 0.05,
    Kp_fall: float | None = None,
    tau_s_fall: float | None = None,
    zeta_fall: float | None = None,
    theta_p_fall: float | None = None,
) -> np.ndarray:
    """
    Causal SOPDT filter on a sampled series (dead time then 2nd-order lag).

    Used as an outlet/instrument overlay (COMSOL ym_* GE), not as PDE inlet forcing.

    Optional ``*_fall`` parameters enable asymmetric valve-close dynamics: rising
    uses delayed command ``y(t-θ_rise)``; falling uses ``y(t-θ_fall)`` (or immediate
    ``y(t)`` if ``theta_p_fall`` is 0) with falling (Kp, τ, ζ). Upper clip against
    the target applies only while rising so purge decay is not snapped to a step.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if t.size == 0:
        return y.copy()
    out = np.zeros_like(y)
    y_state = 0.0
    dy = 0.0
    theta_r = max(float(theta_p), 0.0)
    use_fall = any(v is not None for v in (Kp_fall, tau_s_fall, zeta_fall, theta_p_fall))
    Kp_f = float(Kp if Kp_fall is None else Kp_fall)
    tau_f = float(tau_s if tau_s_fall is None else tau_s_fall)
    zeta_f = float(zeta if zeta_fall is None else zeta_fall)
    theta_f = max(float(theta_p if theta_p_fall is None else theta_p_fall), 0.0)
    t0 = float(t[0])
    for i in range(len(t)):
        t_i = float(t[i])
        dt = float(t[i] - t[i - 1]) if i else 0.0
        u_imm = float(y[i])
        u_del_r = _interp_u(t, y, t_i - theta_r, t0)
        rising = u_imm > u_del_r + 1e-15
        if rising or not use_fall:
            Kp_i, tau_i, zeta_i = float(Kp), float(tau_s), float(zeta)
            u = u_del_r
        else:
            Kp_i, tau_i, zeta_i = Kp_f, tau_f, zeta_f
            u = u_imm if theta_f <= 0.0 else _interp_u(t, y, t_i - theta_f, t0)
        target = Kp_i * u
        inv_tau2 = 1.0 / max(tau_i * tau_i, 1e-30)
        two_zeta_tau = 2.0 * zeta_i / max(tau_i, 1e-30)
        if dt <= 0.0:
            out[i] = y_state
            continue
        n_sub = max(1, int(np.ceil(dt / max(dt_sub, 1e-6))))
        h = dt / n_sub
        for _ in range(n_sub):
            ddy = inv_tau2 * (target - y_state) - two_zeta_tau * dy
            dy = dy + h * ddy
            y_new = y_state + h * dy
            if rising and target > 0.0 and y_new > target:
                # no overshoot above rising setpoint; do not snap falling edges
                y_state = target
                if dy > 0.0:
                    dy = 0.0
            else:
                y_state = max(y_new, 0.0)
        out[i] = y_state
    return out


def _params_from_block(block: dict, row: dict, direction_key: str) -> dict[str, float] | None:
    """Extract Kp/tau/zeta/theta from rising|falling or flat row."""
    directed = row.get(direction_key) if isinstance(row.get(direction_key), dict) else None
    src = directed
    if src is None and direction_key == "rising":
        if any(k in row for k in ("tau_s", "zeta", "theta_p", "Kp")):
            src = row
        elif not row:
            src = block
        else:
            return None
    if src is None:
        return None
    return {
        "Kp": float(src.get("Kp", block.get("Kp", 1.0))),
        "tau_s": float(src.get("tau_s", block.get("tau_s", 1.0))),
        "zeta": float(src.get("zeta", block.get("zeta", 1.0))),
        "theta_p": float(src.get("theta_p", block.get("theta_p", 0.0))),
    }


def _outlet_species_row(block: dict, species: str) -> dict | None:
    """
    Return per_species row for outlet filter, or None if this species is not filtered.

    Legacy: if ``per_species`` is omitted/empty, CO2 may use top-level block defaults.
    If ``per_species`` is non-empty, only listed species are filtered (no CO2 fallback).
    """
    per = block.get("per_species")
    if per is None:
        per = {}
    if not per:
        return {} if species == "CO2" else None
    row = per.get(species)
    return row if isinstance(row, dict) else None


def outlet_filter_params(cfg: DfmConfig, species: str) -> dict[str, float] | None:
    """Return rising Kp, tau_s, zeta, theta_p for outlet overlay, or None if disabled."""
    z0 = (cfg.gas_phase.get("boundaries_mass") or {}).get("z_equals_0") or {}
    block = z0.get("outlet_sopdt_filter") or {}
    if not block.get("enabled", False):
        return None
    row = _outlet_species_row(block, species)
    if row is None:
        return None
    return _params_from_block(block, row, "rising")


def outlet_filter_params_asymmetric(
    cfg: DfmConfig, species: str
) -> tuple[dict[str, float], dict[str, float] | None] | None:
    """Return (rising, falling|None) outlet SOPDT params, or None if disabled."""
    z0 = (cfg.gas_phase.get("boundaries_mass") or {}).get("z_equals_0") or {}
    block = z0.get("outlet_sopdt_filter") or {}
    if not block.get("enabled", False):
        return None
    row = _outlet_species_row(block, species)
    if row is None:
        return None
    rising = _params_from_block(block, row, "rising")
    if rising is None:
        return None
    falling = _params_from_block(block, row, "falling")
    return rising, falling


def filter_outlet_series(
    cfg: DfmConfig,
    species: str,
    t: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """Apply config outlet_sopdt_filter to series y(t); returns y unchanged if disabled."""
    pack = outlet_filter_params_asymmetric(cfg, species)
    if pack is None:
        return np.asarray(y, dtype=float)
    rising, falling = pack
    kwargs: dict[str, float] = dict(rising)
    if falling is not None:
        kwargs.update(
            {
                "Kp_fall": falling["Kp"],
                "tau_s_fall": falling["tau_s"],
                "zeta_fall": falling["zeta"],
                "theta_p_fall": falling["theta_p"],
            }
        )
    return apply_sopdt_filter(t, y, **kwargs)


def initial_gas_concentrations(cfg: DfmConfig, species: list[str]) -> dict[str, float]:
    sim = cfg.raw.get("simulation") or {}
    ic = sim.get("initial_C_mol_m3")
    if ic:
        return {sp: float(ic[sp]) for sp in species}
    if sim.get("cycle"):
        return {sp: 0.0 for sp in species}
    return CycleScheduler(cfg).inlet_at(0.0)
