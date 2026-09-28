"""
CAPE-OPEN DFM atom-ledger outputs (step 5).

Copies compute_balance_summary fields onto flat output parameters.
Does not recompute R_E from CAPE ports (adsorbed q is not on the
Material Object).
"""

from __future__ import annotations

from typing import Any

from adr.cape.contract import BALANCE_TOL_DEFAULT, CapeContractError

_ELEMENTS = ("C", "H", "O")


def require_dfm_ledger(summary: Any) -> dict[str, Any]:
    if not isinstance(summary, dict):
        raise CapeContractError("run JSON has no DFM balance_summary")
    if str(summary.get("schema") or "").startswith("fermentation"):
        raise CapeContractError("fermentation balance_summary is not a DFM ledger")
    if "atom_residual_mol" not in summary or "atom_residual_relative" not in summary:
        raise CapeContractError("balance_summary is not a DFM atom ledger")
    return summary


def ledger_balance_ok(summary: dict[str, Any], tol: float) -> bool:
    rel = summary.get("atom_residual_relative") or {}
    try:
        vals = [abs(float(rel[el])) for el in _ELEMENTS]
    except (KeyError, TypeError, ValueError):
        return False
    return bool(vals) and max(vals) < float(tol)


def _opt_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def dfm_outputs_from_summary(
    summary: dict[str, Any],
    *,
    tol: float = BALANCE_TOL_DEFAULT,
) -> dict[str, Any]:
    """
    Flat SI outputs for ICapeParameter.

    Values are copied from the engine ledger. Port mole flows are not used.
    """
    summary = require_dfm_ledger(summary)
    rel = summary.get("atom_residual_relative") or {}
    mol = summary.get("atom_residual_mol") or {}
    conv = summary.get("conversion") or {}
    hf = summary.get("holdup_final_mol") or {}
    h0 = summary.get("holdup_initial_mol") or {}
    ok = ledger_balance_ok(summary, tol)
    max_rel = None
    try:
        max_rel = max(abs(float(rel[el])) for el in _ELEMENTS)
    except (KeyError, TypeError, ValueError):
        max_rel = None
    return {
        "R_C_rel": _opt_float(rel.get("C")),
        "R_H_rel": _opt_float(rel.get("H")),
        "R_O_rel": _opt_float(rel.get("O")),
        "R_C_mol": _opt_float(mol.get("C")),
        "R_H_mol": _opt_float(mol.get("H")),
        "R_O_mol": _opt_float(mol.get("O")),
        "Y_CH4": _opt_float(conv.get("Y_CH4_per_CO2_in")),
        "n_CH4": _opt_float(conv.get("n_CH4_produced_mol")),
        "N_CO2_ads": _opt_float(hf.get("CO2_ads")),
        "N_CO2_ads_0": _opt_float(h0.get("CO2_ads")),
        "H2O_over_CH4": _opt_float(conv.get("H2O_over_CH4")),
        "X_H2_to_CH4": _opt_float(conv.get("X_H2_to_CH4")),
        "X_CO2_converted": _opt_float(conv.get("X_CO2_converted")),
        "balance_ok": ok,
        "balance_tol": float(tol),
        "max_R_E_rel": max_rel,
        "note": (
            "Copied from DFM compute_balance_summary (includes adsorbed q). "
            "Not recomputed from CAPE ports."
        ),
    }


def enforce_balance_gate(outputs: dict[str, Any], *, fail_if_unbalanced: bool) -> None:
    if not fail_if_unbalanced:
        return
    if outputs.get("balance_ok"):
        return
    max_rel = outputs.get("max_R_E_rel")
    tol = outputs.get("balance_tol")
    raise CapeContractError(
        f"balance_ok is false: max |R_E|/scale={max_rel} >= tol={tol}"
    )


def force_dfm_balance_enabled(raw: dict[str, Any]) -> None:
    """CAPE Calculate() always keeps the DFM atom ledger on."""
    sim = raw.setdefault("simulation", {})
    bal = sim.get("balance")
    if isinstance(bal, dict):
        bal["enabled"] = True
    else:
        sim["balance"] = {"enabled": True}
