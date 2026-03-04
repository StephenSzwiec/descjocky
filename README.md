# DescJocky

A cheminformatics molecular descriptor calculator with pluggable backends.

DescJocky operates in two phases:

- **Phase 1 (Geometry Optimization):** Takes SMILES strings, generates 3-D conformers via RDKit ETKDG, and optimizes geometries with [xtb](https://github.com/grimme-lab/xtb) — each molecule in its own isolated sub-process.
- **Phase 2 (Descriptor Calculation):** Runs pluggable descriptor backends in parallel over the optimized structures, producing a single CSV.

## Architecture

```
SMILES file
    │
    ▼
┌──────────────────────────┐
│  Phase 1: Geometry Opt   │   ProcessPoolExecutor
│  SMILES → ETKDG → xtb    │   one xtb sub-process per molecule
│  Output: SDF per mol     │   temp dirs prevent filename collisions
└──────────────────────────┘
    │
    │  MolRecord (mol_id, smiles, sdf_path)
    │  SDF files are the interchange format
    ▼
┌──────────────────────────┐
│  Phase 2: Descriptors    │   Per-backend parallelism
│                          │
│  ┌────────┐ ┌─────────┐  │   concurrency_safe=True  → ProcessPool
│  │ RDKit  │ │ Mordred │  │   concurrency_safe=False → sequential
│  └────────┘ └─────────┘  │
│  ┌────────┐ ┌─────────┐  │   Backends load Mol from SDF on worker
│  │ Pybel  │ │ Native  │  │   side — no pickling of Mol objects
│  └────────┘ └─────────┘  │
│  ┌─────────────────────┐ │
│  │ Your backend here   │ │
│  └─────────────────────┘ │
└──────────────────────────┘
    │
    ▼
descriptors.csv
```

### Key Design Decisions

**SDF as interchange format.** RDKit Mol objects don't pickle reliably with 3-D conformers across `spawn`-mode process boundaries. Instead, Phase 1 writes one SDF per molecule, and Phase 2 workers reconstitute the Mol on their side via `MolRecord.load_mol()`.

**Per-backend concurrency strategy.** Pure Python backends (Mordred, RDKit) run in a `ProcessPoolExecutor`. Backends that spawn their own subprocesses (Java-based tools, for example) should declare `concurrency_safe = False` and run sequentially to avoid over-subscription.

**Plugin architecture.** Every backend is a subclass of `Backend` with four methods: `setup()`, `compute()`, `teardown()`, and `available()`. Drop a new file in `descjocky/backends/` and it auto-registers on import.

## Backends

| Backend  | Descriptors | 3-D | Dependencies     | Notes                                    |
|----------|-------------|-----|------------------|------------------------------------------|
| RDKit    | ~210        | No  | *(core)*         | Always available. Built-in RDKit descs.  |
| Native   | ~25+        | Yes | *(core)*         | Clean-room implementation. Growing.      |
| Mordred  | ~1,800      | Yes | `mordred`        | Most comprehensive open-source library.  |
| Pybel    | ~20         | No  | `openbabel`      | OpenBabel's perception model.            |

List installed backends:
```
descjocky --list-backends
```

## Installation

```bash
git clone https://github.com/stephenszwiec/descjocky
cd descjocky
conda env create -f environment.yml
conda activate descjocky
```

Or install just the core (no conda needed, but you'll need xtb on PATH):
```bash
pip install -e .               # RDKit + Native only
pip install -e ".[mordred]"    # + Mordred
pip install -e ".[all]"        # + Mordred + Pybel
```

## Usage

### Config file
```ini
[files]
input_file = mols.txt
mol_dir = ./mols
csv_output = descriptors.csv

[settings]
num_workers = 4
xtb_timeout = 600
remove_temp_files = false
backends = all
```

```bash
descjocky -c config.ini
```

### Command line
```bash
descjocky --smiles mols.txt --mol-dir ./mols --backends RDKit,Mordred --workers 8
```

### Skip Phase 1
If you already have optimized SDFs in `mol_dir/optimized/`:
```bash
descjocky --smiles mols.txt --mol-dir ./mols --skip-phase1
```

### Python API
```python
from descjocky import Pipeline

records = Pipeline({
    "input_file": "mols.txt",
    "mol_dir": "./mols",
    "csv_output": "descriptors.csv",
    "num_workers": 4,
    "backends": "RDKit,Native",
}).run()
```

## Writing a Backend

```python
from descjocky.core.backend import Backend
from descjocky.core.mol_record import MolRecord
from descjocky.core.registry import BackendRegistry

class MyBackend(Backend):
    name = "MyDescs"
    concurrency_safe = True
    supports_3d = True

    def setup(self):
        # Called once per worker process
        pass

    def compute(self, record: MolRecord) -> dict[str, float | str]:
        mol = record.load_mol()  # RDKit Mol with 3-D coords
        return {
            "my_descriptor": some_calculation(mol),
        }

    @classmethod
    def available(cls) -> bool:
        return True

BackendRegistry.register(MyBackend)
```

Save this as `descjocky/backends/my_backend.py` and it will be auto-discovered.

## Project Structure

```
descjocky/
├── pyproject.toml              # Packaging, deps, entry point
├── environment.yml             # Conda env (includes xtb)
├── descjocky/
│   ├── __init__.py             # Public API
│   ├── cli.py                  # Argument parsing, config, entry point
│   ├── core/
│   │   ├── backend.py          # Backend ABC (the plugin contract)
│   │   ├── registry.py         # Auto-discovery and registration
│   │   ├── mol_record.py       # Data transfer object (picklable)
│   │   ├── geometry.py         # Phase 1: conformer gen + xtb
│   │   ├── descriptors.py      # Phase 2: parallel backend dispatch
│   │   ├── pipeline.py         # Top-level orchestrator
│   │   └── writer.py           # CSV output
│   └── backends/
│       ├── rdkit_backend.py    # Always available (~210 descs)
│       ├── native_backend.py   # Clean-room calculator (~25+ descs)
│       ├── mordred_backend.py  # Optional (~1800 descs)
│       └── pybel_backend.py    # Optional (~20 descs)
└── tests/
    └── test_core.py
```

## Roadmap

The Native backend is where DescJocky becomes more than an aggregator. 

Planned descriptor families:

- **Topological:** Kier-Hall chi/kappa, Balaban J, information content indices
- **Electronic:** CPSA (charged partial surface area), DPSA, PEOE
- **3-D geometric:** WHIM, GETAWAY, RDF, 3D-MoRSE, autocorrelation
- **Fingerprint-derived:** bit statistics from Morgan, MACCS, atom-pair FPs
- **Pharmacophoric:** pharmacophore-type counts and distances

## References

- **xtb:** Bannwarth, Caldeweyher, Ehlert et al. *WIREs Comput. Mol. Sci.* 2020, 11, e01493.
- **RDKit:** Open-source cheminformatics. [rdkit.org](https://rdkit.org)
- **Mordred:** Moriwaki et al. *J Cheminform* 10, 4 (2018).
- **OpenBabel:** O'Boyle, Morley, Hutchison. *Chemistry Central Journal* 2, 5 (2008).
- **BACE-1 dataset:** Subramanian et al. *J. Chem. Inf. Model.* 2016, 56(10), 1936–1949.
