"""
DescJocky — a cheminformatics molecular-descriptor calculator.

Public API
----------
Pipeline        High-level orchestrator (Phase 1 → Phase 2 → CSV).
BackendRegistry Discover / register descriptor backends at runtime.
MolRecord       Lightweight data-transfer object that flows between phases.
"""

from descjocky.core.pipeline import Pipeline
from descjocky.core.registry import BackendRegistry
from descjocky.core.mol_record import MolRecord

__all__ = ["Pipeline", "BackendRegistry", "MolRecord"]
__version__ = "0.2.0"
