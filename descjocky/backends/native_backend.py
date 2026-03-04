"""
Skeleton: native DescJocky descriptor backend.

This file is a template for the "dream bigger" goal — a clean-room
descriptor calculator implemented directly against RDKit's Mol object
and NumPy, without depending on any third-party descriptor library.

The idea: implement the most impactful descriptors (constitutional,
topological, electronic, geometric) in pure Python/NumPy so that
DescJocky can produce a rich feature set with *zero* optional
dependencies beyond rdkit + numpy.

This is where DescJocky stops being an aggregator and starts being
a descriptor calculator in its own right.

To activate, remove the ``return False`` from ``available()`` and
implement ``compute()``.
"""

from __future__ import annotations

import logging
from typing import ClassVar, Optional

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdmolops

from descjocky.core.backend import Backend
from descjocky.core.mol_record import MolRecord
from descjocky.core.registry import BackendRegistry

log = logging.getLogger(__name__)


class NativeBackend(Backend):
    """Clean-room descriptor calculator — the future of DescJocky.

    Categories to implement (roughly ordered by value-for-effort):

    1. Constitutional
       Atom/bond counts, molecular weight, hydrogen counts,
       ring counts, rotatable bonds, etc.  Many of these overlap
       with RDKit's built-ins, but having our own implementation
       means we control the perception rules exactly.

    2. Topological
       Wiener index, Balaban J, Zagreb indices, Randic connectivity,
       Kier-Hall chi/kappa indices.  Computable from the molecular
       graph alone (no 3-D needed).

    3. Electronic / Charge
       Gasteiger charges → partial charge statistics (mean, max, min,
       sum of positive, sum of negative, DPSA descriptors).
       Requires only 2-D + Gasteiger.

    4. Geometrical (3-D)
       Principal moments of inertia, asphericity, eccentricity,
       radius of gyration, WHIM, GETAWAY, RDF, 3D-MoRSE.
       These require the optimised geometry from Phase 1.

    5. Fingerprint-derived
       Bit counts and statistics from Morgan, MACCS, topological
       fingerprints.  Useful as ML features directly.
    """

    name: ClassVar[str] = "Native"
    concurrency_safe: ClassVar[bool] = True
    supports_3d: ClassVar[bool] = True

    def setup(self) -> None:
        pass

    def compute(self, record: MolRecord) -> dict[str, float | str]:
        mol = record.load_mol()
        if mol is None:
            raise ValueError(record.error or "could not load mol")

        descs: dict[str, float | str] = {}

        # ── 1. Constitutional ──────────────────────────────────────
        descs["n_atoms_heavy"] = mol.GetNumHeavyAtoms()
        descs["n_atoms_total"] = mol.GetNumAtoms()
        descs["n_bonds"] = mol.GetNumBonds()
        descs["n_rings"] = rdMolDescriptors.CalcNumRings(mol)
        descs["n_aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
        descs["n_rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
        descs["n_hba"] = rdMolDescriptors.CalcNumHBA(mol)
        descs["n_hbd"] = rdMolDescriptors.CalcNumHBD(mol)
        descs["mw"] = Descriptors.ExactMolWt(mol)
        descs["logp_wildman_crippen"] = Descriptors.MolLogP(mol)
        descs["tpsa"] = Descriptors.TPSA(mol)
        descs["fraction_csp3"] = rdMolDescriptors.CalcFractionCSP3(mol)

        # ── 2. Topological (graph-based) ───────────────────────────
        dm = rdmolops.GetDistanceMatrix(mol)
        descs["wiener_index"] = float(np.sum(dm) / 2)
        descs["diameter"] = float(np.max(dm))

        # Zagreb indices
        degrees = np.array([
            atom.GetDegree() for atom in mol.GetAtoms()
        ], dtype=float)
        descs["zagreb_m1"] = float(np.sum(degrees ** 2))
        descs["zagreb_m2"] = float(sum(
            degrees[b.GetBeginAtomIdx()] * degrees[b.GetEndAtomIdx()]
            for b in mol.GetBonds()
        ))

        # ── 3. Charge descriptors ──────────────────────────────────
        try:
            Chem.rdPartialCharges.ComputeGasteigerCharges(mol)
            charges = np.array([
                float(atom.GetDoubleProp("_GasteigerCharge"))
                for atom in mol.GetAtoms()
            ])
            # Replace NaN/inf from Gasteiger failures
            charges = np.where(np.isfinite(charges), charges, 0.0)
            descs["charge_mean"] = float(np.mean(charges))
            descs["charge_max"] = float(np.max(charges))
            descs["charge_min"] = float(np.min(charges))
            descs["charge_sum_pos"] = float(np.sum(charges[charges > 0]))
            descs["charge_sum_neg"] = float(np.sum(charges[charges < 0]))
        except Exception:
            for key in ("charge_mean", "charge_max", "charge_min",
                        "charge_sum_pos", "charge_sum_neg"):
                descs[key] = "NA"

        # ── 4. Geometrical (3-D) ──────────────────────────────────
        conf = mol.GetConformer(0) if mol.GetNumConformers() > 0 else None
        if conf is not None:
            coords = np.array(conf.GetPositions())
            centroid = coords.mean(axis=0)
            centered = coords - centroid

            # Radius of gyration
            descs["radius_of_gyration"] = float(
                np.sqrt(np.mean(np.sum(centered ** 2, axis=1)))
            )

            # Inertia tensor → principal moments
            inertia = np.zeros((3, 3))
            for c in centered:
                inertia += np.eye(3) * np.dot(c, c) - np.outer(c, c)
            eigvals = np.sort(np.linalg.eigvalsh(inertia))
            descs["pmi1"] = float(eigvals[0])
            descs["pmi2"] = float(eigvals[1])
            descs["pmi3"] = float(eigvals[2])

            # Normalized PMI ratios (Sauer & Schwarz style)
            if eigvals[2] > 0:
                descs["npr1"] = float(eigvals[0] / eigvals[2])
                descs["npr2"] = float(eigvals[1] / eigvals[2])
            else:
                descs["npr1"] = "NA"
                descs["npr2"] = "NA"

            # Asphericity
            trace = float(eigvals.sum())
            if trace > 0:
                descs["asphericity"] = float(
                    1.5 * np.sum((eigvals - trace / 3) ** 2) / (trace ** 2)
                )
            else:
                descs["asphericity"] = "NA"
        else:
            for key in ("radius_of_gyration", "pmi1", "pmi2", "pmi3",
                        "npr1", "npr2", "asphericity"):
                descs[key] = "NA"

        return descs

    @classmethod
    def available(cls) -> bool:
        # Flip to True when this backend is ready for use
        return True

    @classmethod
    def descriptor_count_hint(cls) -> Optional[int]:
        return 25  # will grow as we add categories


BackendRegistry.register(NativeBackend)
