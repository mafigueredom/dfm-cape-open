"""
CAPE-OPEN DFM reports (step 6).

ICapeUnitReport payload: atom-ledger text and optional outlet C(t).
Never attached to the Product Material Object. C(z)/q(z) are not exported.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

from adr.cape.contract import SPECIES

REPORT_ATOM = "atom_balance"
REPORT_C_RAW = "outlet_C_raw"
REPORT_C_SOPDT = "outlet_C_sopdt"
REPORT_PROFILES = "profiles_C_z"


def _series_ok(series: Any) -> bool:
    return (
        isinstance(series, dict)
        and isinstance(series.get("times_s"), list)
        and series["times_s"]
        and isinstance(series.get("C_mol_m3"), dict)
        and series["C_mol_m3"]
    )


def format_outlet_csv(series: dict[str, Any]) -> str:
    """CSV text: t_s, C_CO2, C_H2, C_CH4, C_H2O, C_N2 [mol/m3]."""
    if not _series_ok(series):
        return ""
    times = series["times_s"]
    conc = series["C_mol_m3"]
    cols = [sp for sp in SPECIES if sp in conc]
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["t_s", *[f"C_{sp}_mol_m3" for sp in cols]])
    n = len(times)
    for i in range(n):
        row = [times[i]]
        for sp in cols:
            vals = conc[sp]
            row.append(vals[i] if i < len(vals) else "")
        writer.writerow(row)
    return buf.getvalue()


def attach_reports(result: dict[str, Any]) -> dict[str, Any]:
    """Fill ``result['reports']`` from ledger text and C(t) series already on the result."""
    names: list[str] = [REPORT_ATOM]
    catalog: dict[str, Any] = {
        REPORT_ATOM: {
            "title": "DFM atom / mass balance",
            "mime": "text/plain",
            "text": str(result.get("report_text") or ""),
            "note": "Whole-reactor ledger including adsorbed CO2. Not a stream.",
        },
        REPORT_PROFILES: {
            "title": "Spatial profiles C(z), q(z)",
            "available": False,
            "note": "Not stored in the current DFM export. Report-only if added later.",
        },
    }
    if _series_ok(result.get("outlet_C_raw")):
        names.append(REPORT_C_RAW)
        catalog[REPORT_C_RAW] = {
            "title": "Outlet mean C(t) raw",
            "mime": "text/csv",
            "available": True,
            "ref": "outlet_C_raw",
            "n_samples": len(result["outlet_C_raw"]["times_s"]),
            "note": "Report series. Not applied to product F_ss.",
        }
    c_t = str(result.get("C_t_export") or "")
    if c_t in ("sopdt", "both"):
        if _series_ok(result.get("outlet_C_sopdt")):
            names.append(REPORT_C_SOPDT)
            catalog[REPORT_C_SOPDT] = {
                "title": "Outlet mean C(t) SOPDT (instrument)",
                "mime": "text/csv",
                "available": True,
                "ref": "outlet_C_sopdt",
                "n_samples": len(result["outlet_C_sopdt"]["times_s"]),
                "note": "Report series from outlet_sopdt_filter. Not applied to product F_ss.",
            }
        else:
            catalog[REPORT_C_SOPDT] = {
                "title": "Outlet mean C(t) SOPDT (instrument)",
                "available": False,
                "note": result.get("outlet_C_sopdt_note")
                or "SOPDT C(t) not filled (needs adapter filter).",
            }
    result["reports"] = {
        "names": names,
        "catalog": catalog,
        "note": (
            "ICapeUnitReport only. Product Material Object stays raw cycle-average."
        ),
    }
    return result


def get_report(result: dict[str, Any], name: str) -> str:
    """Return report body as text (COM GetReport)."""
    if name == REPORT_ATOM:
        cat = (result.get("reports") or {}).get("catalog") or {}
        return str((cat.get(REPORT_ATOM) or {}).get("text") or result.get("report_text") or "")
    if name == REPORT_C_RAW:
        return format_outlet_csv(result.get("outlet_C_raw") or {})
    if name == REPORT_C_SOPDT:
        return format_outlet_csv(result.get("outlet_C_sopdt") or {})
    if name == REPORT_PROFILES:
        return "Spatial profiles C(z), q(z) are not in the current DFM export.\n"
    raise KeyError(f"unknown CAPE report {name!r}")


def write_cape_reports(
    result: dict[str, Any],
    directory: Path | str,
) -> list[Path]:
    """Write atom_balance.txt and outlet_C_*.csv under ``directory``."""
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    attach_reports(result)
    names = list((result.get("reports") or {}).get("names") or [])
    mapping = {
        REPORT_ATOM: d / "atom_balance.txt",
        REPORT_C_RAW: d / "outlet_C_raw.csv",
        REPORT_C_SOPDT: d / "outlet_C_sopdt.csv",
    }
    for name in names:
        path = mapping.get(name)
        if path is None:
            continue
        body = get_report(result, name)
        if not body:
            continue
        path.write_text(body, encoding="utf-8")
        written.append(path)
        cat = result["reports"]["catalog"].setdefault(name, {})
        cat["file"] = str(path)
    return written
