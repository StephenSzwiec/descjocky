"""
BackendRegistry — discover, register, and instantiate descriptor backends.

Backends are discovered in two ways:

1. **Built-in scanning**: everything in ``descjocky.backends`` that
   subclasses ``Backend`` is auto-registered on import.
2. **Explicit registration**: third-party code can call
   ``BackendRegistry.register(MyBackend)`` at any time.

Only backends whose ``available()`` class-method returns ``True`` (i.e.
whose dependencies are installed) are offered to the pipeline.
"""

from __future__ import annotations

import importlib
import logging
import pkgutil
from typing import Iterator, Type

from descjocky.core.backend import Backend

log = logging.getLogger(__name__)


class BackendRegistry:
    """Singleton-ish registry of descriptor backends."""

    _backends: dict[str, Type[Backend]] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    @classmethod
    def register(cls, backend_cls: Type[Backend]) -> None:
        """Register a backend class (idempotent)."""
        name = backend_cls.name
        if name in cls._backends:
            return
        cls._backends[name] = backend_cls
        log.debug("Registered backend: %s", name)

    @classmethod
    def get(cls, name: str) -> Type[Backend]:
        cls._ensure_discovered()
        return cls._backends[name]

    @classmethod
    def available(cls) -> Iterator[Type[Backend]]:
        """Yield every registered backend whose deps are present."""
        cls._ensure_discovered()
        for backend_cls in cls._backends.values():
            try:
                if backend_cls.available():
                    yield backend_cls
            except Exception:
                log.warning("Backend %s.available() raised; skipping",
                            backend_cls.name, exc_info=True)

    @classmethod
    def names(cls) -> list[str]:
        cls._ensure_discovered()
        return list(cls._backends.keys())

    # ------------------------------------------------------------------
    # Auto-discovery
    # ------------------------------------------------------------------

    _discovered = False

    @classmethod
    def _ensure_discovered(cls) -> None:
        if cls._discovered:
            return
        cls._discovered = True
        # Import every module under descjocky.backends; each module
        # is expected to call register() in its module body.
        package = importlib.import_module("descjocky.backends")
        for finder, module_name, is_pkg in pkgutil.iter_modules(
            package.__path__, prefix="descjocky.backends."
        ):
            try:
                importlib.import_module(module_name)
            except Exception:
                log.debug("Could not import backend module %s",
                          module_name, exc_info=True)
