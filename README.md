# descjocky
##  A Cheminformatics Descriptor Aggregator

DescJocky is a Python-based, parallelized tool for calculating and aggregating molecular descriptors from a variety of open-source cheminformatics libraries.

It operates in two distinct phases:
- **Phase 1 (Geometry Optimization):** Takes a text file of SMILES strings, generates 3D conformers, and runs a single, multi-threaded xtb job to produce high-quality, optimized geometries for each molecule.
- **Phase 2 (Descriptor Calculation):** Uses a multiprocessing pool to calculate descriptors in parallel for each of the optimized structures.

This tool was developed as a non-GUI replacement for the currently inaccessible (ChemDes)[http://www.scbdd.com/chemdes/] platform, with the hopes it is useful to the research community.

Target Libraries:
- PaDEL (via padelpy)
- BlueDesc (via BlueDesc-pywrapper)
- Mordred
- Pybel (OpenBabel)
- Chemopy (via chemopy2)

### Installation
DescJocky is designed to run in a self-contained conda environment due to dependencies to `xtb`.

Clone or download the repository: 
```
git clone https://github.com/stephenszwiec/descjocky
cd descjocky
```

(Or, just download `descjocky.py`, `descjocky.yml`, and `descjocky.ini` into the same directory).

**Create the Conda Environment:** This command reads the `descjocky.yml` file and installs all required dependencies into a new environment named descjocky:
```
conda env create -f descjocky.yml
```
Then, activate this environment any time you want to run the script.
```
conda activate descjocky
```

**Configuration:** You can then set up your own config .ini file as follows:
```
[files]
# REQUIRED: Path to your input SMILES file
input_file = /path/to/your/input.txt
# REQUIRED: A directory for intermediate files
mol_dir = /path/to/directory 
# OPTIONAL: Path for the final CSV. Defaults to 'descriptors.csv' in the current directory.
csv_output_file = /path/to/descriptors.csv

[settings]
# OPTIONAL: Max parallel workers. 0 = use all cores.
num_workers = 0
# OPTIONAL: Set to true to delete all intermediate 3D files.
remove_temp_files = false
```

**Usage:** 
While in your descjocky environment, run the python script.
```
python descjocky.py -c descjocky.ini
```
Or, if a file `config.ini` is in the same directory: `python descjocky.py`

**Results:** 
The script will log its progress to the console (stderr) and create the descriptors.csv (or your specified output file) when complete.

Note: column descriptor output is expected to be overlapping and redundant. Preprocess your CSV to remove highly colinear/correlated columns as needed.

### References 

- **xtb:** Semiempirical Extended Tight-Binding `xtb` lives (on GitHub)[https://github.com/grimme-lab/xtb]
  - C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme WIREs Comput. Mol. Sci., 2020, 11, e01493. DOI: 10.1002/wcms.1493
- **rdkit:** Open-source cheminformatics (https://rdkit.org)[https://rdkit.org] also (on GitHub)[https://github.com/rdkit]
- **mordred:** a molecular descriptor calculator (on GitHub)[https://github.com/mordred-descriptor/mordred] 
  - Moriwaki, H., Tian, YS., Kawashita, N. et al. Mordred: a molecular descriptor calculator. J Cheminform 10, 4 (2018). https://doi.org/10.1186/s13321-018-0258-y
- **PaDELPy:** A python wrapper for PaDEL-Descriptor (on GitHub)[https://github.com/ecrl/padelpy]
  - Yap CW. PaDEL-descriptor: an open source software to calculate molecular descriptors and fingerprints. J Comput Chem. 2011 May;32(7):1466-74. doi: 10.1002/jcc.21707. Epub 2010 Dec 17. PMID: 21425294.
- **BlueDesc-pywrapper:** A python wrapper for BlueDesc molecular descriptors (on GitHub)[https://github.com/OlivierBeq/BlueDesc_pywrapper]
  - BlueDesc Descriptor Calculator (last update: 03-10-08) originally by Georg Hinselmann and released under  GNU GPL
- **Pybel:** a Pythonic API for `openbabel` (on GitHub)[https://github.com/openbabel/openbabel]
  - O'Boyle, N.M., Morley, C. & Hutchison, G.R. Pybel: a Python wrapper for the OpenBabel cheminformatics toolkit. Chemistry Central Journal 2, 5 (2008). https://doi.org/10.1186/1752-153X-2-5
- **chemopy2** a Python library calcuating molecular descriptors (on GitHub)[https://github.com/OlivierBeq/chemopy/]
  - Cao et al., ChemoPy: freely available python package for computational biology and chemoinformatics. Bioinformatics 2013; 29(8), 1092–1094. doi:10.1093/bioinformatics/btt105
- **BACE-1** The chemical data included in this repository was taken from the following source:
  - Subramanian, G.; Ramsundar, B.; Pande, V.; Denny, R.A. Computational Modeling of β-Secretase 1 (BACE-1) Inhibitors Using Ligand Based Approaches. J. Chem. Inf. Model. 2016, 56 (10), 1936–1949. DOI: 10.1021/acs.jcim.6b00290
