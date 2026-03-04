"""
Pybel (OpenBabel) backend — ~20 descriptors from a C++ cheminformatics
library with a different perception model than RDKit.

Including Pybel gives users access to descriptors computed with
OpenBabel's SMARTS/aromaticity model, which can be complementary
to RDKit's, especially for ML feature sets where diversity of
representation helps.

Note: Pybel works from SMILES (2-D only).  The optimised geometry
is not used here, but the backend could be extended to read the SDF
for 3-D-aware OB descriptors in the future.
"""

from __future__ import annotations

import logging
from typing import ClassVar, Optional

from descjocky.core.backend import Backend
from descjocky.core.mol_record import MolRecord
from descjocky.core.registry import BackendRegistry

log = logging.getLogger(__name__)


class PybelBackend(Backend):
    name: ClassVar[str] = "Pybel"
    concurrency_safe: ClassVar[bool] = True
    supports_3d: ClassVar[bool] = False

    def setup(self) -> None:
        # Importing pybel can be slow (loads the C++ lib); do it once.
        from openbabel import pybel  # noqa: F401

    def compute(self, record: MolRecord) -> dict[str, float | str]:
        from openbabel import pybel

        mol = pybel.readstring("smi", record.smiles)
        descs = mol.calcdesc()
        # calcdesc() returns a dict of str→float already
        return {k: v for k, v in descs.items()}

    @classmethod
    def available(cls) -> bool:
        try:
            from openbabel import pybel  # noqa: F401
            return True
        except ImportError:
            return False

    @classmethod
    def descriptor_count_hint(cls) -> Optional[int]:
        return 20


BackendRegistry.register(PybelBackend)
