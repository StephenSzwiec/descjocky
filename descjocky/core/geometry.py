"""
Phase 1 — 3-D conformer generation and geometry optimisation with xtb.

Each molecule gets:
    1. SMILES → RDKit embedded 3-D conformer (ETKDG v3).
    2. Written to an individual SDF in ``initial/``.
    3. Optimised via ``xtb --opt`` in its own temporary working directory,
       so output filenames (``xtbopt.sdf``) never collide.
    4. The optimised SDF is moved to ``optimized/<mol_id>.sdf``.

Concurrency
-----------
Geometry optimisation is CPU-bound *and* subprocess-bound.  We use
``concurrent.futures.ProcessPoolExecutor`` with each task running
xtb in its own working directory.  The xtb ``-P`` flag is set to 1
so each xtb process is single-threaded; parallelism comes from
running many of them simultaneously.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Sequence

from rdkit import Chem
from rdkit.Chem import AllChem

from descjocky.core.mol_record import MolRecord

log = logging.getLogger(__name__)

# xtb threads per subprocess — always 1; we parallelise at the
# process level instead.
_XTB_THREADS = 1


def _embed_conformer(smiles: str) -> Chem.Mol | None:
    """SMILES → RDKit Mol with a 3-D conformer (ETKDG v3)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        # Fallback: try without the distance-geometry seed constraint
        status = AllChem.EmbedMolecule(mol)
    if status != 0:
        return None
    # Quick MMFF pre-optimisation so xtb has a reasonable starting point
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass  # non-fatal; xtb will handle it
    return mol


def _write_initial_sdf(mol: Chem.Mol, path: Path) -> None:
    """Write a single Mol to an SDF file, preserving hydrogens."""
    with Chem.SDWriter(str(path)) as w:
        w.write(mol)


def _run_xtb(
    sdf_in: Path,
    sdf_out: Path,
    xtb_path: str,
    timeout: int | None = 600,
) -> bool:
    """Run ``xtb --opt`` on *sdf_in*, writing the result to *sdf_out*.

    Each invocation gets its own temporary working directory so the
    hard-coded ``xtbopt.sdf`` output does not collide.

    Returns True on success.
    """
    with tempfile.TemporaryDirectory(prefix="xtb_") as tmpdir:
        # xtb reads from the input path but writes into cwd
        cmd = [
            xtb_path,
            str(sdf_in),
            "--opt",
            "-P", str(_XTB_THREADS),
        ]
        try:
            result = subprocess.run(
                cmd,
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            log.warning("xtb timed out for %s", sdf_in.name)
            return False

        opt_sdf = Path(tmpdir) / "xtbopt.sdf"
        if result.returncode != 0 or not opt_sdf.exists():
            log.debug("xtb stderr for %s:\n%s", sdf_in.name,
                      result.stderr[-500:] if result.stderr else "(empty)")
            return False

        shutil.move(str(opt_sdf), str(sdf_out))
        return True


def _optimise_one(
    record: MolRecord,
    initial_dir: Path,
    opt_dir: Path,
    xtb_path: str,
    timeout: int | None,
) -> MolRecord:
    """End-to-end: embed → write SDF → xtb → move result."""
    mol = _embed_conformer(record.smiles)
    if mol is None:
        record.error = "conformer embedding failed"
        return record

    mol.SetProp("_Name", f"{record.mol_id}: {record.smiles}")
    init_sdf = initial_dir / f"{record.mol_id}.sdf"
    _write_initial_sdf(mol, init_sdf)

    opt_sdf = opt_dir / f"{record.mol_id}.sdf"
    ok = _run_xtb(init_sdf, opt_sdf, xtb_path, timeout=timeout)
    if ok:
        record.sdf_path = opt_sdf
    else:
        record.error = "xtb optimisation failed"
    return record


def run_geometry_optimisation(
    records: Sequence[MolRecord],
    *,
    mol_dir: Path,
    xtb_path: str,
    num_workers: int = 4,
    timeout: int | None = 600,
) -> list[MolRecord]:
    """Phase 1: parallel geometry optimisation of all molecules.

    Parameters
    ----------
    records
        MolRecords populated with SMILES only.
    mol_dir
        Root directory for intermediate files.
    xtb_path
        Resolved path to the ``xtb`` executable.
    num_workers
        Max parallel xtb subprocesses.
    timeout
        Per-molecule xtb wall-clock timeout in seconds.

    Returns
    -------
    list[MolRecord]
        The same records, now with ``sdf_path`` set on success
        or ``error`` set on failure.
    """
    initial_dir = mol_dir / "initial"
    opt_dir = mol_dir / "optimized"
    initial_dir.mkdir(parents=True, exist_ok=True)
    opt_dir.mkdir(parents=True, exist_ok=True)

    log.info(
        "Phase 1: optimising %d molecules with %d workers",
        len(records), num_workers,
    )

    results: list[MolRecord] = []

    with ProcessPoolExecutor(max_workers=num_workers) as pool:
        futures = {
            pool.submit(
                _optimise_one, rec, initial_dir, opt_dir, xtb_path, timeout
            ): rec.mol_id
            for rec in records
        }
        for future in as_completed(futures):
            mol_id = futures[future]
            try:
                result = future.result()
            except Exception as exc:
                # Should not happen — _optimise_one catches internally —
                # but be defensive.
                log.error("Unexpected error for %s: %s", mol_id, exc)
                # Find the original record and mark it
                for rec in records:
                    if rec.mol_id == mol_id:
                        rec.error = str(exc)
                        results.append(rec)
                        break
                continue
            if result.ok:
                log.info("Optimised %s", mol_id)
            else:
                log.warning("Failed %s: %s", mol_id, result.error)
            results.append(result)

    success = sum(1 for r in results if r.ok)
    log.info("Phase 1 complete: %d / %d succeeded", success, len(records))
    return results
