"""
Backend ABC — the contract every descriptor backend must fulfil.

Design goals
------------
* Thin interface: just ``name``, ``setup``, ``compute``, ``available``.
* Backends are free to accept *either* an RDKit Mol or an SDF path;
  the ``MolRecord`` carries both.
* Backends return a plain ``dict[str, float|str]``; the pipeline
  handles prefixing and merging.
* Backends that shell out to Java / Fortran / etc. can override
  ``concurrency_safe`` to signal they need subprocess isolation
  rather than in-process multiprocessing.
"""

from __future__ import annotations

import abc
from typing import ClassVar, Optional

from descjocky.core.mol_record import MolRecord


class Backend(abc.ABC):
    """Abstract base for a descriptor-calculation backend."""

    # Human-readable name; used as the column-prefix in the output CSV.
    name: ClassVar[str]

    # If False, the pipeline will run this backend in isolated
    # subprocesses rather than a shared multiprocessing pool.
    concurrency_safe: ClassVar[bool] = True

    # If True, the backend can make use of 3-D coordinates.
    supports_3d: ClassVar[bool] = False

    # ---- lifecycle ----------------------------------------------------

    def setup(self) -> None:
        """One-time initialisation (heavy object creation, JVM start, …).

        Called once per *worker process*, not once globally.
        """

    def teardown(self) -> None:
        """Optional cleanup (temp dirs, JVM shutdown, …)."""

    # ---- core ---------------------------------------------------------

    @abc.abstractmethod
    def compute(self, record: MolRecord) -> dict[str, float | str]:
        """Calculate descriptors for a single molecule.

        Parameters
        ----------
        record : MolRecord
            Carries the SMILES, the optimized-SDF path, and a helper
            ``load_mol()`` to get an RDKit Mol with 3-D coords.

        Returns
        -------
        dict
            Descriptor-name → value.  Keys should be *unprefixed*;
            the pipeline adds the backend ``name`` as a prefix.

        Raises
        ------
        Exception
            Any exception is caught by the pipeline, logged, and the
            molecule is marked with an error rather than crashing the run.
        """
        ...

    # ---- introspection ------------------------------------------------

    @classmethod
    def available(cls) -> bool:
        """Return True if this backend's dependencies are installed.

        The default implementation tries to import the backend module
        itself; override for more nuanced checks (e.g. Java on PATH).
        """
        return True

    @classmethod
    def descriptor_count_hint(cls) -> Optional[int]:
        """Approximate number of descriptors this backend produces.

        Used for progress-bar sizing and nothing else.
        """
        return None
