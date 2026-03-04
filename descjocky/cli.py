"""
CLI entry point for DescJocky.

Usage::

    descjocky -c config.ini
    descjocky --smiles mols.txt --mol-dir ./mols --backends RDKit,Mordred
    descjocky --list-backends
"""

from __future__ import annotations

import argparse
import configparser
import logging
import multiprocessing as mp
import sys
from pathlib import Path

from descjocky.core.pipeline import Pipeline
from descjocky.core.registry import BackendRegistry


def _setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] (%(processName)s) %(message)s",
        handlers=[
            logging.StreamHandler(sys.stderr),
            logging.FileHandler("descjocky.log", mode="w"),
        ],
    )


def _parse_ini(path: str) -> dict:
    """Parse a .ini config file into a flat dict for the Pipeline."""
    cp = configparser.ConfigParser()
    cp.read(path)
    config: dict = {}

    if "files" in cp:
        config["input_file"] = cp["files"].get("input_file", "mols.txt")
        config["mol_dir"] = cp["files"].get("mol_dir", "./mols")
        config["csv_output"] = cp["files"].get("csv_output", "descriptors.csv")

    if "settings" in cp:
        raw_workers = int(cp["settings"].get("num_workers", "0"))
        config["num_workers"] = raw_workers if raw_workers > 0 else mp.cpu_count()
        config["remove_temp"] = cp["settings"].getboolean("remove_temp_files", False)
        config["xtb_timeout"] = int(cp["settings"].get("xtb_timeout", "600"))
        config["backends"] = cp["settings"].get("backends", "all")

    return config


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="descjocky",
        description="DescJocky — a cheminformatics molecular-descriptor calculator",
    )
    p.add_argument("-c", "--config", default=None,
                   help="Path to .ini config file")
    p.add_argument("--smiles", default=None,
                   help="Path to input SMILES file (overrides config)")
    p.add_argument("--mol-dir", default=None,
                   help="Working directory for SDF files (overrides config)")
    p.add_argument("--csv-output", default=None,
                   help="Output CSV path (overrides config)")
    p.add_argument("--backends", default=None,
                   help="Comma-separated backend names, or 'all'")
    p.add_argument("--workers", type=int, default=None,
                   help="Number of parallel workers")
    p.add_argument("--skip-phase1", action="store_true",
                   help="Skip geometry optimisation (use existing SDFs)")
    p.add_argument("--list-backends", action="store_true",
                   help="List available backends and exit")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    _setup_logging(verbose=args.verbose)

    # Handle --list-backends
    if args.list_backends:
        print("Available descriptor backends:\n")
        for backend_cls in BackendRegistry.available():
            hint = backend_cls.descriptor_count_hint()
            n = f"~{hint} descriptors" if hint else "unknown count"
            safe = "concurrent" if backend_cls.concurrency_safe else "sequential"
            threed = "2D+3D" if backend_cls.supports_3d else "2D only"
            print(f"  {backend_cls.name:12s}  {n:20s}  {threed:8s}  ({safe})")
        print()
        sys.exit(0)

    # Build config: start from INI, then overlay CLI args
    if args.config:
        config = _parse_ini(args.config)
    elif Path("config.ini").exists():
        config = _parse_ini("config.ini")
    else:
        config = {}

    # CLI overrides
    if args.smiles:
        config["input_file"] = args.smiles
    if args.mol_dir:
        config["mol_dir"] = args.mol_dir
    if args.csv_output:
        config["csv_output"] = args.csv_output
    if args.backends:
        config["backends"] = args.backends
    if args.workers is not None:
        config["num_workers"] = args.workers
    if args.skip_phase1:
        config["skip_phase1"] = True

    # Validate we have the minimum required config
    if "input_file" not in config:
        parser.error("No input SMILES file specified. Use -c config.ini or --smiles path")
    config.setdefault("mol_dir", "./mols")

    print("DescJocky — Molecular Descriptor Calculator")
    print("=" * 45)

    pipeline = Pipeline(config)
    pipeline.run()


if __name__ == "__main__":
    try:
        mp.set_start_method("spawn")
    except RuntimeError:
        pass
    main()
