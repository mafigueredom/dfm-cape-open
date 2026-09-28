"""Export / report view over FieldStateHistory."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from adr.model.state import FieldStateHistory
from fenicsx.balance_diagnostics import format_balance_report
from fenicsx.run_export import save_reactor_state_json


class ResultExportView:
    def export_json(
        self,
        state: FieldStateHistory,
        *,
        path: Path | str,
        species: list[str],
        config_path: Path | str | None = None,
        cfg: Any = None,
        meta: dict[str, Any] | None = None,
    ) -> Path:
        if meta is None and cfg is not None:
            raw = cfg.raw.get("meta") or {}
            meta = {
                "fenics_implementation": raw.get("fenics_implementation"),
                "alignment_target": raw.get("alignment_target"),
                "factory_id": getattr(state, "factory_id", None)
                or raw.get("factory_id"),
            }
        return save_reactor_state_json(
            path,
            state,
            species=species,
            config_path=config_path,
            meta=meta,
        )

    def print_balance(self, state: FieldStateHistory) -> None:
        if state.balance_summary is None:
            return
        schema = str(state.balance_summary.get("schema") or "")
        if schema.startswith("fermentation"):
            print("  Fermentation summary:", state.balance_summary.get("note", ""))
            return
        print(format_balance_report(state.balance_summary))
