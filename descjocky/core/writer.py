"""
Output writer — deterministic CSV with consistent column ordering.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Sequence

from descjocky.core.mol_record import MolRecord

log = logging.getLogger(__name__)


def write_csv(records: Sequence[MolRecord], out_path: Path) -> int:
    """Write descriptor results to a CSV.

    Parameters
    ----------
    records
        All MolRecords (including failed ones — they are included
        with an ``error`` column for traceability).
    out_path
        Destination CSV path.

    Returns
    -------
    int
        Number of rows written.
    """
    if not records:
        log.error("No records to write.")
        return 0

    # Build ordered column list:  mol_id, smiles, error, then
    # all descriptor columns in sorted order.
    desc_keys: set[str] = set()
    for rec in records:
        desc_keys.update(rec.descriptors.keys())

    columns = ["mol_id", "smiles", "error"] + sorted(desc_keys)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, restval="NA",
                                extrasaction="ignore")
        writer.writeheader()
        for rec in records:
            row = {
                "mol_id": rec.mol_id,
                "smiles": rec.smiles,
                "error": rec.error or "",
                **rec.descriptors,
            }
            writer.writerow(row)

    log.info("Wrote %d rows (%d columns) to %s",
             len(records), len(columns), out_path)
    return len(records)
