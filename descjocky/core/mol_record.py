"""
MolRecord — lightweight, picklable data-transfer object.

Every molecule in the pipeline is represented as a MolRecord.  The key
design choice: the *SDF file path* is the canonical representation that
crosses process boundaries.  RDKit Mol objects are reconstituted lazily
on the worker side, avoiding pickle issues entirely.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Optional

from rdkit import Chem


@dataclasses.dataclass(slots=True)
class MolRecord:
    """One molecule moving through the descriptor pipeline.

    Attributes
    ----------
    mol_id : str
        Unique, human-readable identifier (e.g. ``mol_0001``).
    smiles : str
        Original input SMILES.
    sdf_path : Path | None
        Path to the optimized SDF file after Phase 1.
        ``None`` if geometry optimization has not yet run or failed.
    descriptors : dict[str, float | str]
        Accumulated descriptor values.  Populated during Phase 2.
    error : str | None
        If a phase fails for this molecule, a human-readable reason.
    """

    mol_id: str
    smiles: str
    sdf_path: Optional[Path] = None
    descriptors: dict = dataclasses.field(default_factory=dict)
    error: Optional[str] = None

    # ------------------------------------------------------------------
    # Lazy Mol access — only on the worker that needs it
    # ------------------------------------------------------------------

    def load_mol(self, remove_hs: bool = False) -> Optional[Chem.Mol]:
        """Read the optimized SDF and return an RDKit Mol (with 3-D coords).

        Returns ``None`` and sets ``self.error`` if loading fails.
        """
        if self.sdf_path is None or not self.sdf_path.exists():
            self.error = f"SDF not found: {self.sdf_path}"
            return None
        try:
            supplier = Chem.SDMolSupplier(str(self.sdf_path), removeHs=remove_hs)
            mol = supplier[0]
            if mol is None:
                self.error = f"RDKit could not parse {self.sdf_path}"
            return mol
        except Exception as exc:
            self.error = f"SDF load error: {exc}"
            return None

    @property
    def ok(self) -> bool:
        """True when no error has been recorded."""
        return self.error is None
