"""CLI: cape_request.json → cape_result.json (no COM)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_V1 = Path(__file__).resolve().parents[2]
if str(_V1) not in sys.path:
    sys.path.insert(0, str(_V1))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run the DFM CAPE adapter (cape_request_v1 → cape_result_v1)."
    )
    parser.add_argument("request", type=Path, help="cape_request_v1 JSON")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Write cape_result_v1 JSON (default: stdout)",
    )
    parser.add_argument("--t-end", type=float, default=None, help="Override t_end_s")
    parser.add_argument("--dt", type=float, default=None, help="Override dt_s")
    parser.add_argument(
        "--export-run",
        type=Path,
        default=None,
        help="Also write the fenics_run JSON",
    )
    parser.add_argument("--summary", action="store_true", help="Print solver summary")
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=None,
        help="Write atom_balance.txt and outlet_C_*.csv (ICapeUnitReport files)",
    )
    args = parser.parse_args(argv)

    from adr.cape.adapter import run_cape_request
    from adr.cape.contract import dump_cape_json, load_cape_request

    req = load_cape_request(args.request)
    result = run_cape_request(
        req,
        base_dir=_V1,
        t_end=args.t_end,
        dt=args.dt,
        export_json=args.export_run,
        print_summary=args.summary,
        reports_dir=args.reports_dir,
    )
    if args.output is not None:
        dump_cape_json(result, args.output)
    else:
        json.dump(result, sys.stdout, indent=2, ensure_ascii=False)
        sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
