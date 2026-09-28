"""Abstract Factory: one concrete factory per physical ADR system."""

from __future__ import annotations

from adr.factories.base import ProblemFactory
from adr.factories.registry import get_factory, list_factories, register_factory

__all__ = [
    "ProblemFactory",
    "get_factory",
    "list_factories",
    "register_factory",
]
