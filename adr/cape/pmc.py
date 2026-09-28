"""
CAPE-OPEN 1.1 Unit Operation wrapper (step 9), COM-free Calculate().

The Windows/Wine COM PMC (``com/DfmMethanation``) copies SI numbers from
Material Objects onto this mapping, runs Docker, and copies the result back.
This module is the tested source of that mapping. It does not implement
ICapeUnit in-process and does not recompute ``R_C`` from ports.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

from adr.cape.contract import (
    ALLOWED_FACTORY_ID,
    CAS_TO_SPECIES,
    REQUEST_SCHEMA,
    SPECIES,
    SPECIES_TO_CAS,
    CapeContractError,
    parse_cape_request,
)

# COM identity (must match com/DfmMethanation/Guids.cs).
CLSID = "{187B8CDA-40ED-463E-B455-4D45A7D209FD}"
PROGID = "PhD.DFM.Methanation.1"
UNIT_NAME = "DFM Methanation Reactor"
CAPE_VERSION = "1.1"
CATID_UNIT = "{678C09A5-7D66-11D2-A67D-00105A42887F}"
CATID_PMC = "{678C09A1-7D66-11D2-A67D-00105A42887F}"
CATID_CONSUMES_THERMO = "{4150C28A-EE06-403F-A871-87AFEC38A249}"
CATID_THERMO_10 = "{0D562DC8-EA8E-4210-AB39-B66513C0CD09}"
CATID_THERMO_11 = "{4667023A-5A8E-4CCA-AB6D-9D78C5112FED}"

PORT_FEEDS = (
    ("Feed_ads", "ads", True),
    ("Feed_purge", "purge", True),
    ("Feed_rxn", "rxn", True),
    ("Feed_purge2", "purge2", False),
)
PORT_PRODUCT = "Product"
OUTPUT_PARAM_NAMES = (
    "R_C_rel",
    "R_H_rel",
    "R_O_rel",
    "Y_CH4",
    "n_CH4",
    "N_CO2_ads",
    "balance_ok",
)
DOCKER_IMAGE_DEFAULT = "dfm-methanation-cape:v1"
_Y_MAPPED_MIN = 1e-12


def _as_float_list(value: Any) -> list[float]:
    if value is None:
        return []
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, str):
        return [float(value)]
    return [float(x) for x in value]


def _as_str_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(x) for x in value]


def y_from_property_package(
    cas_list: list[str],
    mole_fractions: list[float],
    *,
    feed_id: str,
) -> tuple[dict[str, float], list[str]]:
    """Map PP compounds by CAS. Unmapped keys are dropped and mapped y is renormalized."""
    notes: list[str] = []
    y = {sp: 0.0 for sp in SPECIES}
    unmapped: list[str] = []
    n = min(len(cas_list), len(mole_fractions))
    for i in range(n):
        cas = str(cas_list[i]).strip()
        frac = float(mole_fractions[i])
        sp = CAS_TO_SPECIES.get(cas) or (cas if cas in SPECIES else None)
        if sp is None:
            if frac > _Y_MAPPED_MIN:
                unmapped.append(cas)
            continue
        y[sp] += frac
    if unmapped:
        notes.append(
            f"feeds.{feed_id}: ignored unmapped PP compounds {unmapped} (not created)"
        )
    total = sum(y.values())
    if total <= _Y_MAPPED_MIN:
        raise CapeContractError(
            f"feeds.{feed_id}: no mapped CAPE-OPEN compounds "
            f"(need CAS {list(CAS_TO_SPECIES)})"
        )
    if abs(total - 1.0) > 1e-6:
        y = {sp: val / total for sp, val in y.items()}
        notes.append(f"feeds.{feed_id}: renormalized mapped y (sum was {total:g})")
    return y, notes


def require_product_slate(cas_list: list[str]) -> list[str]:
    """PP must own all five DFM compounds so the product MO can be written."""
    have = {str(c).strip() for c in cas_list}
    missing = [cas for cas in SPECIES_TO_CAS.values() if cas not in have]
    if missing:
        raise CapeContractError(
            f"property package missing CAS {missing}; cannot write Product"
        )
    return list(cas_list)


def material_to_feed(
    material: dict[str, Any] | None,
    feed_id: str,
    *,
    required: bool,
) -> tuple[dict[str, Any] | None, list[str]]:
    """
    SI snapshot of a connected Material Object → convention-B feed.

    Reads ``F_mol_s``, ``P_Pa``, ``y`` (CAS or species). Stream ``T_K`` is
    recorded on the snapshot but is **not** used as ``T_isothermal_K``.
    """
    if material is None:
        if required:
            raise CapeContractError(f"port for feeds.{feed_id} is not connected")
        return None, []
    cas = _as_str_list(material.get("cas") or material.get("CAS"))
    y_raw = material.get("y")
    notes: list[str] = []
    if isinstance(y_raw, dict):
        keys = [str(k) for k in y_raw]
        y, n1 = y_from_property_package(
            keys, [float(y_raw[k]) for k in y_raw], feed_id=feed_id
        )
        notes.extend(n1)
    elif cas and material.get("mole_fractions") is not None:
        y, n1 = y_from_property_package(
            cas, _as_float_list(material.get("mole_fractions")), feed_id=feed_id
        )
        notes.extend(n1)
    else:
        raise CapeContractError(f"feeds.{feed_id} needs y or (cas + mole_fractions)")
    flows = _as_float_list(material.get("F_mol_s") or material.get("total_flow_mol_s"))
    if not flows:
        raise CapeContractError(f"feeds.{feed_id}: molar flow is missing")
    press = _as_float_list(material.get("P_Pa") or material.get("pressure_Pa"))
    if not press or press[0] <= 0.0:
        raise CapeContractError(f"feeds.{feed_id}: pressure must be > 0 Pa")
    t_stream = material.get("T_K") or material.get("temperature_K")
    feed = {
        "id": feed_id,
        "F_mol_s": float(flows[0]),
        "P_Pa": float(press[0]),
        "y": y,
    }
    if t_stream is not None:
        feed["T_stream_K"] = float(_as_float_list(t_stream)[0])
        notes.append(
            f"feeds.{feed_id}: stream T ignored for T_bed "
            "(energy_thermal.T_isothermal_K)"
        )
    return feed, notes


def request_from_socket(socket: dict[str, Any]) -> dict[str, Any]:
    """Build ``cape_request_v1`` from a COM socket snapshot (no FEniCS)."""
    notes: list[str] = []
    materials = socket.get("feeds") or socket.get("materials") or {}
    feeds: dict[str, Any] = {}
    required = {fid: req for _name, fid, req in PORT_FEEDS}
    for _name, fid, req in PORT_FEEDS:
        feed, n1 = material_to_feed(materials.get(fid), fid, required=req)
        notes.extend(n1)
        feeds[fid] = feed
    t_purge2 = (
        ((socket.get("parameters") or {}).get("cycle") or {}).get("t_purge2")
    )
    if t_purge2 is not None and float(t_purge2) <= 0.0:
        feeds["purge2"] = None
    product_cas = _as_str_list(
        socket.get("product_cas") or socket.get("pp_cas") or []
    )
    if product_cas:
        require_product_slate(product_cas)
    raw = {
        "schema": REQUEST_SCHEMA,
        "factory_id": socket.get("factory_id") or ALLOWED_FACTORY_ID,
        "config_path": socket.get("config_path")
        or "config/dfm_config.pow_v07_10_60_opt-1.json",
        "C_t_export": socket.get("C_t_export") or "both",
        "balance_tol": socket.get("balance_tol", 0.001),
        "fail_if_unbalanced": bool(socket.get("fail_if_unbalanced", False)),
        "feeds": feeds,
        "parameters": socket.get("parameters") or {},
    }
    req = parse_cape_request(raw)
    req["pmc_notes"] = notes
    req["product_cas"] = product_cas
    return req


def product_material_from_result(
    result: dict[str, Any],
    *,
    pp_cas: list[str] | None = None,
) -> dict[str, Any]:
    """Cycle-average raw product for the Material Object. SOPDT is never applied."""
    product = result.get("product") or {}
    if "outlet_C_sopdt" in product:
        raise CapeContractError("SOPDT series must not be on the product stream")
    f_ss = product.get("F_ss_mol_s") or {}
    y = product.get("y") or {}
    total = sum(float(f_ss.get(sp, 0.0)) for sp in SPECIES)
    cas_list = list(pp_cas or [SPECIES_TO_CAS[sp] for sp in SPECIES])
    require_product_slate(cas_list)
    y_cas: dict[str, float] = {}
    f_cas: dict[str, float] = {}
    for cas in cas_list:
        sp = CAS_TO_SPECIES.get(cas)
        if sp is None:
            y_cas[cas] = 0.0
            f_cas[cas] = 0.0
            continue
        y_cas[cas] = float(y.get(sp, 0.0))
        f_cas[cas] = float(f_ss.get(sp, 0.0))
    return {
        "T_K": float(product.get("T_K") or 593.15),
        "P_Pa": float(product.get("P_Pa") or 101325.0),
        "F_mol_s": total,
        "y": y_cas,
        "F_ss_mol_s": f_cas,
        "flash": "TP",
        "note": product.get("note") or "SOPDT not applied",
    }


def output_parameters_from_result(result: dict[str, Any]) -> dict[str, Any]:
    """Copy engine ledger outputs. Do not recompute from ports."""
    out = dict(result.get("outputs") or {})
    return {name: out.get(name) for name in OUTPUT_PARAM_NAMES}


def report_names_from_result(result: dict[str, Any]) -> list[str]:
    names = list(((result.get("reports") or {}).get("names")) or [])
    return names


def run_docker_engine(
    request: dict[str, Any],
    work_dir: Path | str,
    *,
    image: str | None = None,
    extra_args: list[str] | None = None,
    docker_bin: str = "docker",
) -> dict[str, Any]:
    """Write request JSON, ``docker run`` the step-8 image, return cape_result_v1."""
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    req_path = work / "cape_request.json"
    res_path = work / "cape_result.json"
    payload = {k: v for k, v in request.items() if k not in ("pmc_notes", "product_cas")}
    req_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    img = image or os.environ.get("DFM_CAPE_IMAGE") or DOCKER_IMAGE_DEFAULT
    uid = os.getuid() if hasattr(os, "getuid") else 0
    gid = os.getgid() if hasattr(os, "getgid") else 0
    cmd = [
        docker_bin,
        "run",
        "--rm",
        "-u",
        f"{uid}:{gid}",
        "-v",
        f"{work.resolve()}:/data",
        img,
        "/data/cape_request.json",
        "-o",
        "/data/cape_result.json",
        *(extra_args or []),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise CapeContractError(
            "Docker CAPE engine failed "
            f"(exit {proc.returncode}): {proc.stderr or proc.stdout}"
        )
    if not res_path.is_file():
        raise CapeContractError("engine did not write cape_result.json")
    return json.loads(res_path.read_text(encoding="utf-8"))


def run_local_engine(
    request: dict[str, Any],
    *,
    base_dir: Path | str | None = None,
    t_end: float | None = None,
    dt: float | None = None,
) -> dict[str, Any]:
    """In-process FEniCS path (Linux tests). Not used by the COM DLL."""
    from adr.cape.adapter import run_cape_request

    return run_cape_request(
        request,
        base_dir=base_dir,
        t_end=t_end,
        dt=dt,
        print_summary=False,
    )


def calculate(
    socket: dict[str, Any],
    *,
    engine: str = "docker",
    work_dir: Path | str | None = None,
    result: dict[str, Any] | None = None,
    base_dir: Path | str | None = None,
    t_end: float | None = None,
    dt: float | None = None,
    extra_args: list[str] | None = None,
) -> dict[str, Any]:
    """
    ``ICapeUnit.Calculate()`` without COM.

    ``engine``: ``docker`` (COM default), ``local`` (dolfinx), or ``result``
    (apply a precomputed cape_result; no solve).
    """
    req = request_from_socket(socket)
    if engine == "result":
        if result is None:
            raise CapeContractError("engine=result requires a cape_result_v1 object")
        cap = result
    elif engine == "local":
        cap = run_local_engine(req, base_dir=base_dir, t_end=t_end, dt=dt)
    elif engine == "docker":
        work = Path(work_dir or (Path.cwd() / "cape_pmc_work"))
        args = list(extra_args or [])
        if t_end is not None:
            args.extend(["--t-end", str(t_end)])
        if dt is not None:
            args.extend(["--dt", str(dt)])
        cap = run_docker_engine(req, work, extra_args=args or None)
    else:
        raise CapeContractError(f"unknown CAPE engine {engine!r}")

    pp_cas = req.get("product_cas") or _as_str_list(socket.get("product_cas"))
    if not pp_cas:
        pp_cas = [SPECIES_TO_CAS[sp] for sp in SPECIES]
    product = product_material_from_result(cap, pp_cas=pp_cas)
    outputs = output_parameters_from_result(cap)
    names = report_names_from_result(cap)
    return {
        "request": req,
        "result": cap,
        "product": product,
        "outputs": outputs,
        "reports": names,
        "notes": list(req.get("pmc_notes") or []),
    }


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="COM-free ICapeUnit.Calculate() (socket JSON → cape_result)."
    )
    parser.add_argument("socket", type=Path, help="cape_socket_v1 or cape_request_v1 JSON")
    parser.add_argument("-o", "--output", type=Path, default=None)
    parser.add_argument(
        "--engine",
        choices=("docker", "local", "result"),
        default=os.environ.get("CAPE_PMC_ENGINE", "docker"),
    )
    parser.add_argument("--result", type=Path, default=None, help="Precomputed cape_result")
    parser.add_argument("--work-dir", type=Path, default=None)
    parser.add_argument("--t-end", type=float, default=None)
    parser.add_argument("--dt", type=float, default=None)
    args = parser.parse_args(argv)

    raw = json.loads(args.socket.read_text(encoding="utf-8"))
    if raw.get("schema") == REQUEST_SCHEMA:
        socket = {
            "feeds": raw.get("feeds") or {},
            "parameters": raw.get("parameters") or {},
            "config_path": raw.get("config_path"),
            "C_t_export": raw.get("C_t_export"),
            "balance_tol": raw.get("balance_tol", 0.001),
            "fail_if_unbalanced": raw.get("fail_if_unbalanced", False),
            "factory_id": raw.get("factory_id"),
        }
        for fid, feed in list(socket["feeds"].items()):
            if isinstance(feed, dict) and "y" in feed:
                socket["feeds"][fid] = feed
    else:
        socket = raw
    pre = None
    if args.result is not None:
        pre = json.loads(args.result.read_text(encoding="utf-8"))
        args.engine = "result"
    out = calculate(
        socket,
        engine=args.engine,
        work_dir=args.work_dir,
        result=pre,
        t_end=args.t_end,
        dt=args.dt,
    )
    text = json.dumps(out, indent=2, ensure_ascii=False)
    if args.output is not None:
        args.output.write_text(text + "\n", encoding="utf-8")
    else:
        sys.stdout.write(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
