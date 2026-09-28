"""Registry of ProblemFactory implementations by factory_id."""

from __future__ import annotations

from typing import Any, Type

from adr.factories.base import ProblemFactory

_REGISTRY: dict[str, Type[ProblemFactory]] = {}


def register_factory(cls: Type[ProblemFactory]) -> Type[ProblemFactory]:
    fid = str(getattr(cls, "factory_id", "") or "").strip()
    if not fid:
        raise ValueError(f"{cls.__name__} missing factory_id")
    _REGISTRY[fid] = cls
    return cls


def _ensure_builtins() -> None:
    if _REGISTRY:
        return
    from adr.factories.methanation_dfm import MethanationDfmFactory

    register_factory(MethanationDfmFactory)
    try:
        from adr.factories.fermentation import FermentationFactory

        register_factory(FermentationFactory)
    except ImportError:
        pass
    try:
        from adr.factories.system_c import SystemCFactory

        register_factory(SystemCFactory)
    except ImportError:
        pass


def list_factories() -> list[str]:
    _ensure_builtins()
    return sorted(_REGISTRY)


def get_factory(
    factory_id: str | None,
    config: Any,
    *,
    default: str = "methanation_dfm",
) -> ProblemFactory:
    """
    Resolve factory_id from argument or config meta / simulation block.

    Lookup order when factory_id is None:
      simulation.factory_id → meta.factory_id → default
    """
    _ensure_builtins()
    fid = (factory_id or "").strip()
    if not fid:
        raw = getattr(config, "raw", None) or {}
        sim = raw.get("simulation") or {}
        meta = raw.get("meta") or {}
        fid = str(sim.get("factory_id") or meta.get("factory_id") or default).strip()
    if fid not in _REGISTRY:
        known = ", ".join(list_factories())
        raise KeyError(f"Unknown factory_id={fid!r}. Known: {known}")
    return _REGISTRY[fid](config)
