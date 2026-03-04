"""
Pipeline — the top-level orchestrator.

    Pipeline(config).run()

is the entire application.  Config is a plain dict (parsed from INI,
TOML, CLI args, or constructed programmatically).
"""

from __future__ import annotations

import logging
import shutil
import sys
from pathlib import Path
from typing import Sequence, Type

from descjocky.core.backend import Backend
from descjocky.core.descriptors import run_descriptor_calculation
from descjocky.core.geometry import run_geometry_optimisation
from descjocky.core.mol_record import MolRecord
from descjocky.core.registry import BackendRegistry
from descjocky.core.writer import write_csv

log = logging.getLogger(__name__)


class Pipeline:
    """End-to-end descriptor-calculation pipeline.

    Parameters
    ----------
    config : dict
        Expected keys::

            input_file       Path to a text file of SMILES (one per line).
            mol_dir          Working directory for intermediate SDF files.
            csv_output       Output CSV path  [default: descriptors.csv].
            num_workers      Parallel workers  [default: 4].
            xtb_timeout      Per-molecule xtb timeout in seconds [default: 600].
            remove_temp      Delete mol_dir after the run  [default: False].
            backends         List of backend names, or "all"  [default: "all"].
            skip_phase1      If True, skip geometry opt (assume SDFs exist).
    """

    def __init__(self, config: dict) -> None:
        self.config = config
        self._validate()

    def _validate(self) -> None:
        c = self.config
        c.setdefault("csv_output", "descriptors.csv")
        c.setdefault("num_workers", 4)
        c.setdefault("xtb_timeout", 600)
        c.setdefault("remove_temp", False)
        c.setdefault("backends", "all")
        c.setdefault("skip_phase1", False)

        self.input_file = Path(c["input_file"]).resolve()
        self.mol_dir = Path(c["mol_dir"]).resolve()
        self.csv_output = Path(c["csv_output"]).resolve()
        self.num_workers = int(c["num_workers"])
        self.xtb_timeout = int(c["xtb_timeout"]) if c["xtb_timeout"] else None
        self.remove_temp = bool(c["remove_temp"])
        self.skip_phase1 = bool(c["skip_phase1"])

        if not self.input_file.exists():
            log.fatal("Input file not found: %s", self.input_file)
            sys.exit(1)

        # Resolve xtb
        self.xtb_path = shutil.which("xtb")
        if self.xtb_path is None and not self.skip_phase1:
            log.fatal("xtb not found on PATH")
            sys.exit(1)

    def _resolve_backends(self) -> list[Type[Backend]]:
        """Turn the config 'backends' value into concrete classes."""
        requested = self.config["backends"]

        if requested == "all":
            result = list(BackendRegistry.available())
        else:
            if isinstance(requested, str):
                requested = [s.strip() for s in requested.split(",")]
            result = []
            for name in requested:
                try:
                    cls = BackendRegistry.get(name)
                    if cls.available():
                        result.append(cls)
                    else:
                        log.warning("Backend '%s' requested but not available", name)
                except KeyError:
                    log.warning("Unknown backend '%s'", name)

        if not result:
            log.fatal("No descriptor backends available.  Install at least one "
                      "(e.g. `pip install mordred`).")
            sys.exit(1)

        log.info("Active backends: %s",
                 ", ".join(b.name for b in result))
        return result

    def _read_smiles(self) -> list[MolRecord]:
        """Parse the input file into MolRecords."""
        records = []
        with open(self.input_file) as fh:
            for i, line in enumerate(fh, 1):
                smi = line.strip()
                if not smi or smi.startswith("#"):
                    continue
                mol_id = f"mol_{i:04d}"
                records.append(MolRecord(mol_id=mol_id, smiles=smi))
        if not records:
            log.fatal("Input file is empty: %s", self.input_file)
            sys.exit(1)
        log.info("Read %d SMILES from %s", len(records), self.input_file)
        return records

    def run(self) -> list[MolRecord]:
        """Execute the full pipeline and return the final records."""
        backends = self._resolve_backends()
        records = self._read_smiles()

        # ── Phase 1 ──
        if self.skip_phase1:
            log.info("Skipping Phase 1 (geometry optimisation)")
            opt_dir = self.mol_dir / "optimized"
            for rec in records:
                candidate = opt_dir / f"{rec.mol_id}.sdf"
                if candidate.exists():
                    rec.sdf_path = candidate
                else:
                    rec.error = "pre-optimised SDF not found"
        else:
            records = run_geometry_optimisation(
                records,
                mol_dir=self.mol_dir,
                xtb_path=self.xtb_path,
                num_workers=self.num_workers,
                timeout=self.xtb_timeout,
            )

        # ── Phase 2 ──
        records = run_descriptor_calculation(
            records,
            backends,
            num_workers=self.num_workers,
        )

        # ── Output ──
        write_csv(records, self.csv_output)

        # ── Cleanup ──
        if self.remove_temp:
            log.info("Removing temporary directory: %s", self.mol_dir)
            try:
                shutil.rmtree(self.mol_dir)
            except Exception as exc:
                log.error("Could not remove temp dir: %s", exc)

        ok = sum(1 for r in records if r.ok and r.descriptors)
        log.info("Pipeline complete. %d / %d molecules have descriptors.",
                 ok, len(records))
        return records
