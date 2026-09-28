"""MVC View: plots, JSON export, balance reporting."""

from __future__ import annotations

from adr.view.export import ResultExportView

__all__ = ["ResultExportView", "ResultPlotView"]


def __getattr__(name: str):
    if name == "ResultPlotView":
        from adr.view.plots import ResultPlotView

        return ResultPlotView
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
