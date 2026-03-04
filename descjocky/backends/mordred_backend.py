"""
Mordred backend — ~1,800 molecular descriptors including 3-D.

Mordred is a pure-Python descriptor calculator that works directly
with RDKit Mol objects.  It is concurrency-safe and the single most
feature-rich open-source descriptor library available.
"""

from __future__ import annotations

import logging
from typing import ClassVar, Optional

from descjocky.core.backend import Backend
from descjocky.core.mol_record import MolRecord
from descjocky.core.registry import BackendRegistry

log = logging.getLogger(__name__)


class MordredBackend(Backend):
    name: ClassVar[str] = "Mordred"
    concurrency_safe: ClassVar[bool] = True
    supports_3d: ClassVar[bool] = True

    def __init__(self) -> None:
        self._calc = None

    def setup(self) -> None:
        from mordred import Calculator, descriptors
        self._calc = Calculator(descriptors, ignore_3D=False)
        log.debug("Mordred calculator initialised with %d descriptors",
                  len(self._calc.descriptors))

    def compute(self, record: MolRecord) -> dict[str, float | str]:
        mol = record.load_mol()
        if mol is None:
            raise ValueError(record.error or "could not load mol")

        result = self._calc(mol)
        # Mordred returns a Result object; convert to dict,
        # coercing errors to NaN strings.
        out: dict[str, float | str] = {}
        for desc, val in zip(self._calc.descriptors, result):
            key = str(desc)
            if isinstance(val, (int, float)):
                out[key] = val
            elif hasattr(val, "error"):
                out[key] = "NA"
            else:
                try:
                    out[key] = float(val)
                except (TypeError, ValueError):
                    out[key] = str(val)
        return out

    @classmethod
    def available(cls) -> bool:
        try:
            import mordred  # noqa: F401
            return True
        except ImportError:
            return False

    @classmethod
    def descriptor_count_hint(cls) -> Optional[int]:
        return 1826


BackendRegistry.register(MordredBackend)
