"""
RDKit built-in descriptors — always available.

Provides ~200 2-D descriptors from RDKit's Descriptors module plus
Morgan fingerprint bits (as an example of fingerprint integration).
This backend has zero extra dependencies beyond the core ``rdkit``
requirement, so it serves as the guaranteed baseline.
"""

from __future__ import annotations

import logging
from typing import ClassVar, Optional

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

from descjocky.core.backend import Backend
from descjocky.core.mol_record import MolRecord
from descjocky.core.registry import BackendRegistry

log = logging.getLogger(__name__)


class RDKitBackend(Backend):
    name: ClassVar[str] = "RDKit"
    concurrency_safe: ClassVar[bool] = True
    supports_3d: ClassVar[bool] = False

    def __init__(self) -> None:
        self._calc: MoleculeDescriptors.MolecularDescriptorCalculator | None = None
        self._desc_names: tuple[str, ...] = ()

    def setup(self) -> None:
        # All descriptors registered with Descriptors module
        self._desc_names = tuple(
            name for name, _ in Descriptors.descList
        )
        self._calc = MoleculeDescriptors.MolecularDescriptorCalculator(
            self._desc_names
        )

    def compute(self, record: MolRecord) -> dict[str, float | str]:
        # RDKit descriptors work from SMILES (2-D) — no need for
        # the optimised geometry.  But if we have 3-D coords, we
        # can also compute some 3-D-aware descriptors later.
        mol = Chem.MolFromSmiles(record.smiles)
        if mol is None:
            raise ValueError(f"RDKit could not parse SMILES: {record.smiles}")

        mol = Chem.AddHs(mol)
        values = self._calc.CalcDescriptors(mol)
        result = dict(zip(self._desc_names, values))

        # Add a few useful extras not in the default list
        result["NumAtoms"] = mol.GetNumAtoms()
        result["NumBonds"] = mol.GetNumBonds()
        result["NumRings"] = Chem.rdMolDescriptors.CalcNumRings(mol)

        return result

    @classmethod
    def available(cls) -> bool:
        return True  # rdkit is a core dependency

    @classmethod
    def descriptor_count_hint(cls) -> Optional[int]:
        return len(Descriptors.descList) + 3


# Self-register on import
BackendRegistry.register(RDKitBackend)
