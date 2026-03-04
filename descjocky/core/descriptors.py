"""
Phase 2 — parallel descriptor calculation across pluggable backends.

Architecture
------------
For each backend, we choose a concurrency strategy:

* ``concurrency_safe = True`` (default): the backend runs in a
  ``ProcessPoolExecutor`` where each worker calls ``setup()`` once,
  then processes a stream of MolRecords.

* ``concurrency_safe = False``: the backend spawns its own
  subprocesses (e.g. Java for PaDEL).  We run it sequentially
  or with limited parallelism to avoid over-subscription.

Results are merged into each MolRecord's ``.descriptors`` dict
with the backend name as a key prefix.
"""

from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Sequence, Type

from descjocky.core.backend import Backend
from descjocky.core.mol_record import MolRecord

log = logging.getLogger(__name__)


# ── Worker helpers (top-level functions for pickling) ──────────────────

# Per-process backend instance, initialised once via _init_worker.
_worker_backend: Backend | None = None


def _init_worker(backend_cls: Type[Backend]) -> None:
    """Initialise the backend once in each worker process."""
    global _worker_backend
    _worker_backend = backend_cls()
    _worker_backend.setup()


def _compute_one(record: MolRecord) -> tuple[str, dict[str, float | str] | None, str | None]:
    """Run the process-local backend on one MolRecord.

    Returns (mol_id, descriptors_or_None, error_or_None).
    """
    assert _worker_backend is not None
    try:
        descs = _worker_backend.compute(record)
        return (record.mol_id, descs, None)
    except Exception as exc:
        return (record.mol_id, None, f"{_worker_backend.name}: {exc}")


# ── Public API ────────────────────────────────────────────────────────

def run_descriptor_calculation(
    records: Sequence[MolRecord],
    backends: Sequence[Type[Backend]],
    *,
    num_workers: int = 4,
) -> list[MolRecord]:
    """Phase 2: calculate descriptors for every (successful) molecule.

    Parameters
    ----------
    records
        MolRecords from Phase 1 (only those with ``record.ok`` are
        processed; the rest are passed through unchanged).
    backends
        Backend *classes* to use (already filtered for availability).
    num_workers
        Max parallel worker processes per concurrency-safe backend.

    Returns
    -------
    list[MolRecord]
        The same list, with ``.descriptors`` populated.
    """
    # Index records by mol_id for fast lookup
    by_id: dict[str, MolRecord] = {r.mol_id: r for r in records}
    ok_records = [r for r in records if r.ok]

    if not ok_records:
        log.warning("No successfully optimised molecules to compute descriptors for.")
        return list(records)

    for backend_cls in backends:
        _run_one_backend(ok_records, backend_cls, by_id, num_workers)

    return list(records)


def _run_one_backend(
    ok_records: list[MolRecord],
    backend_cls: Type[Backend],
    by_id: dict[str, MolRecord],
    num_workers: int,
) -> None:
    """Run a single backend across all molecules."""
    name = backend_cls.name
    log.info("Running backend '%s' on %d molecules …", name, len(ok_records))

    if backend_cls.concurrency_safe:
        _run_concurrent(ok_records, backend_cls, by_id, num_workers)
    else:
        _run_sequential(ok_records, backend_cls, by_id)


def _run_concurrent(
    ok_records: list[MolRecord],
    backend_cls: Type[Backend],
    by_id: dict[str, MolRecord],
    num_workers: int,
) -> None:
    """Run a concurrency-safe backend in a process pool."""
    name = backend_cls.name
    successes = 0

    with ProcessPoolExecutor(
        max_workers=num_workers,
        initializer=_init_worker,
        initargs=(backend_cls,),
    ) as pool:
        futures = {
            pool.submit(_compute_one, rec): rec.mol_id
            for rec in ok_records
        }
        for future in as_completed(futures):
            mol_id = futures[future]
            try:
                mid, descs, err = future.result()
            except Exception as exc:
                log.error("[%s] Unexpected worker error for %s: %s",
                          name, mol_id, exc)
                continue

            rec = by_id[mid]
            if err:
                log.warning("[%s] %s", name, err)
            elif descs:
                # Prefix keys and merge
                rec.descriptors.update(
                    {f"{name}_{k}": v for k, v in descs.items()}
                )
                successes += 1

    log.info("[%s] computed for %d / %d molecules", name, successes, len(ok_records))


def _run_sequential(
    ok_records: list[MolRecord],
    backend_cls: Type[Backend],
    by_id: dict[str, MolRecord],
) -> None:
    """Run a non-concurrency-safe backend one molecule at a time."""
    name = backend_cls.name
    backend = backend_cls()
    backend.setup()
    successes = 0

    try:
        for rec in ok_records:
            try:
                descs = backend.compute(rec)
                if descs:
                    rec.descriptors.update(
                        {f"{name}_{k}": v for k, v in descs.items()}
                    )
                    successes += 1
            except Exception as exc:
                log.warning("[%s] %s: %s", name, rec.mol_id, exc)
    finally:
        backend.teardown()

    log.info("[%s] computed for %d / %d molecules", name, successes, len(ok_records))
