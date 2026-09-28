"""
Whole-reactor molar holdups, atom closure, and conversion metrics.

Axisymmetric volume measure matches the diluted-species weak form weight ``r dx``
plus the physical factor ``2π``:

    ∫_Ω3D φ dV = 2π ∫ φ r dr dz
"""

from __future__ import annotations

import math
from typing import Any

from mpi4py import MPI
from dolfinx import fem
import ufl

from dfm_config import DfmConfig


def _trap_add(cum: float, y0: float, y1: float, dt: float) -> float:
    if dt <= 0.0:
        return cum
    return cum + 0.5 * (y0 + y1) * dt


class HoldupAssembler:
    """Cached scalar forms for gas and adsorbed CO₂ holdups [mol]."""

    def __init__(
        self,
        domain,
        r_coord,
        *,
        epsilon: float,
        rho_a: float,
        comm: MPI.Intracomm | None = None,
    ):
        self._domain = domain
        self._r = r_coord
        self._eps = float(epsilon)
        self._rho_a = float(rho_a)
        self._comm = comm or domain.comm
        self._two_pi = 2.0 * math.pi

    def gas_holdup_mol(self, C: fem.Function) -> float:
        """N_gas = ε ∫ C dV [mol] over the full 3D bed (void inventory)."""
        form = fem.form(self._two_pi * self._eps * C * self._r * ufl.dx)
        local = fem.assemble_scalar(form)
        return float(self._comm.allreduce(local, op=MPI.SUM))

    def ads_holdup_mol(self, q: fem.Function) -> float:
        """N_ads = ρ_a ∫ q dV [mol] (CO₂ on solid, bed volume)."""
        form = fem.form(self._two_pi * self._rho_a * q * self._r * ufl.dx)
        local = fem.assemble_scalar(form)
        return float(self._comm.allreduce(local, op=MPI.SUM))


def sample_holdups(
    assembler: HoldupAssembler,
    C_n: dict[str, fem.Function],
    q_n: fem.Function,
    species: list[str],
) -> dict[str, float]:
    out = {sp: assembler.gas_holdup_mol(C_n[sp]) for sp in species}
    out["CO2_ads"] = assembler.ads_holdup_mol(q_n)
    return out


def atom_inventory(holdup: dict[str, float]) -> dict[str, float]:
    """Element moles in gas + adsorbed CO₂."""
    co2 = float(holdup.get("CO2", 0.0))
    h2 = float(holdup.get("H2", 0.0))
    ch4 = float(holdup.get("CH4", 0.0))
    h2o = float(holdup.get("H2O", 0.0))
    ads = float(holdup.get("CO2_ads", 0.0))
    return {
        "C": co2 + ch4 + ads,
        "H": 2.0 * h2 + 4.0 * ch4 + 2.0 * h2o,
        "O": 2.0 * co2 + h2o + 2.0 * ads,
    }


def atom_feed(F: dict[str, float]) -> dict[str, float]:
    """Element moles from species molar amounts (in or out cumulative)."""
    co2 = float(F.get("CO2", 0.0))
    h2 = float(F.get("H2", 0.0))
    ch4 = float(F.get("CH4", 0.0))
    h2o = float(F.get("H2O", 0.0))
    return {
        "C": co2 + ch4,
        "H": 2.0 * h2 + 4.0 * ch4 + 2.0 * h2o,
        "O": 2.0 * co2 + h2o,
    }


def compute_balance_summary(
    *,
    species: list[str],
    holdup_0: dict[str, float],
    holdup_f: dict[str, float],
    F_in: dict[str, float],
    F_out: dict[str, float],
) -> dict[str, Any]:
    """
    End-of-run closure and conversion.

    Atom residual R_E = ΔN_E − (F_E^in − F_E^out). Ideal ≈ 0.
    Relative residual uses max(|F_in|, |N0|, |Nf|, 1e-30) as scale.
    """
    dN = {sp: float(holdup_f.get(sp, 0.0)) - float(holdup_0.get(sp, 0.0)) for sp in species}
    dN["CO2_ads"] = float(holdup_f.get("CO2_ads", 0.0)) - float(holdup_0.get("CO2_ads", 0.0))

    species_residual = {}
    for sp in species:
        ads = dN["CO2_ads"] if sp == "CO2" else 0.0
        inv = dN[sp] + ads
        flux = float(F_in.get(sp, 0.0)) - float(F_out.get(sp, 0.0))
        # For reactive species this is not expected to be ~0 (stoich exchange).
        species_residual[sp] = inv - flux

    atoms0 = atom_inventory(holdup_0)
    atomsf = atom_inventory(holdup_f)
    Fin_a = atom_feed(F_in)
    Fout_a = atom_feed(F_out)
    atom_residual = {}
    atom_residual_rel = {}
    for el in ("C", "H", "O"):
        d_inv = atomsf[el] - atoms0[el]
        flux = Fin_a[el] - Fout_a[el]
        R = d_inv - flux
        scale = max(abs(Fin_a[el]), abs(atoms0[el]), abs(atomsf[el]), abs(Fout_a[el]), 1e-30)
        atom_residual[el] = R
        atom_residual_rel[el] = R / scale

    n_ch4 = float(F_out.get("CH4", 0.0)) - float(F_in.get("CH4", 0.0))
    n_h2o = float(F_out.get("H2O", 0.0)) - float(F_in.get("H2O", 0.0))
    F_co2_in = float(F_in.get("CO2", 0.0))
    F_h2_in = float(F_in.get("H2", 0.0))
    dH2 = float(F_in.get("H2", 0.0)) - float(F_out.get("H2", 0.0)) - dN.get("H2", 0.0)
    # CO2 consumed from feed that is not left in gas/ads at end nor exited as CO2:
    # carbon to CH4 path: n_CH4 vs CO2 accounted
    co2_accounted = (
        float(F_out.get("CO2", 0.0))
        + float(holdup_f.get("CO2", 0.0))
        + float(holdup_f.get("CO2_ads", 0.0))
        - float(holdup_0.get("CO2", 0.0))
        - float(holdup_0.get("CO2_ads", 0.0))
    )
    co2_to_products = F_co2_in - co2_accounted  # should ≈ n_ch4 if C closes

    conversion = {
        "n_CH4_produced_mol": n_ch4,
        "n_H2O_produced_mol": n_h2o,
        "H2O_over_CH4": (n_h2o / n_ch4) if abs(n_ch4) > 1e-30 else None,
        "Y_CH4_per_CO2_in": (n_ch4 / F_co2_in) if F_co2_in > 1e-30 else None,
        "X_H2_to_CH4": (4.0 * n_ch4 / F_h2_in) if F_h2_in > 1e-30 else None,
        "X_CO2_to_CH4": (n_ch4 / F_co2_in) if F_co2_in > 1e-30 else None,
        "note_X_CO2_to_CH4": "n_CH4 / F_CO2_in (yield on CO2 fed); see also co2_converted_mol",
        "co2_converted_mol": co2_to_products,
        "X_CO2_converted": (co2_to_products / F_co2_in) if F_co2_in > 1e-30 else None,
        "four_CH4_vs_H2_consumed_mol": {
            "4_n_CH4": 4.0 * n_ch4,
            "H2_consumed_incl_holdup": dH2,
        },
    }

    return {
        "holdup_initial_mol": {k: float(v) for k, v in holdup_0.items()},
        "holdup_final_mol": {k: float(v) for k, v in holdup_f.items()},
        "delta_holdup_mol": {k: float(v) for k, v in dN.items()},
        "F_in_mol": {sp: float(F_in.get(sp, 0.0)) for sp in species},
        "F_out_mol": {sp: float(F_out.get(sp, 0.0)) for sp in species},
        "species_residual_mol": species_residual,
        "species_residual_note": (
            "Δ(N_gas[+N_ads for CO2]) − (F_in−F_out). Non-zero for reactive "
            "species is expected; use atom_residual for closure."
        ),
        "atom_inventory_initial_mol": atoms0,
        "atom_inventory_final_mol": atomsf,
        "atom_F_in_mol": Fin_a,
        "atom_F_out_mol": Fout_a,
        "atom_residual_mol": atom_residual,
        "atom_residual_relative": atom_residual_rel,
        "conversion": conversion,
    }


def accumulate_flux_step(
    F: dict[str, float],
    *,
    species: list[str],
    Q0: float,
    Q1: float,
    C0: dict[str, float],
    C1: dict[str, float],
    dt: float,
) -> None:
    """In-place trapezoidal update of cumulative molar flow ∫ Q C dt."""
    for sp in species:
        y0 = float(Q0) * float(C0.get(sp, 0.0))
        y1 = float(Q1) * float(C1.get(sp, 0.0))
        F[sp] = _trap_add(float(F.get(sp, 0.0)), y0, y1, dt)


def format_balance_report(summary: dict[str, Any]) -> str:
    lines = ["--- Mass / atom balance summary ---"]
    Fin = summary["F_in_mol"]
    Fout = summary["F_out_mol"]
    lines.append(
        "  F_in  [mol]: "
        + ", ".join(f"{sp}={Fin[sp]:.6g}" for sp in Fin)
    )
    lines.append(
        "  F_out [mol]: "
        + ", ".join(f"{sp}={Fout[sp]:.6g}" for sp in Fout)
    )
    h0 = summary["holdup_initial_mol"]
    hf = summary["holdup_final_mol"]
    lines.append(
        f"  N_CO2_ads: {h0.get('CO2_ads', 0):.6g} → {hf.get('CO2_ads', 0):.6g} mol"
    )
    ar = summary["atom_residual_mol"]
    arr = summary["atom_residual_relative"]
    lines.append(
        "  Atom residual R=ΔN−(Fin−Fout) [mol] (rel): "
        + ", ".join(f"{el}={ar[el]:.3e} ({arr[el]:.3e})" for el in ("C", "H", "O"))
    )
    conv = summary["conversion"]
    lines.append(
        f"  n_CH4={conv['n_CH4_produced_mol']:.6g} mol, "
        f"H2O/CH4={conv['H2O_over_CH4']}, "
        f"Y_CH4={conv['Y_CH4_per_CO2_in']}, "
        f"X_H2→CH4={conv['X_H2_to_CH4']}, "
        f"X_CO2_converted={conv['X_CO2_converted']}"
    )
    four = conv["four_CH4_vs_H2_consumed_mol"]
    lines.append(
        f"  4·n_CH4={four['4_n_CH4']:.6g}, "
        f"H2_consumed(incl. holdup)={four['H2_consumed_incl_holdup']:.6g} mol"
    )
    return "\n".join(lines)


def balance_enabled(cfg: DfmConfig) -> bool:
    sim = cfg.raw.get("simulation") or {}
    bal = sim.get("balance")
    if bal is None:
        return True
    if isinstance(bal, bool):
        return bal
    return bool(bal.get("enabled", True))


def balance_every_n_steps(cfg: DfmConfig) -> int:
    sim = cfg.raw.get("simulation") or {}
    bal = sim.get("balance") or {}
    if isinstance(bal, dict):
        n = int(bal.get("every_n_steps", 1))
        return max(1, n)
    return 1
