#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DescJocky: A Python-based Cheminformatics Descriptor Aggregator

This script wraps multiple descriptor libraries to calculate and aggregate
molecular descriptors for a list of input SMILES strings, parallelizing
the computation and outputting a single CSV file.
"""


# --- Import Standard Libraries --- 

import os
import sys
import logging
import argparse
import configparser
import csv
import subprocess
import shutil
from pathlib import Path
from functools import partial
import multiprocessing as mp

# --- Import Cheminformatics Libraries ---

from rdkit import Chem
from rdkit.Chem import AllChem
import padelpy
from BlueDesc_pywrapper import BlueDesc
from mordred import Calculator, descriptors
from chemopy import ChemoPy
from chemopy import Fingerprint, Fingerprint3D 
from openbabel import pybel

# --- Global Configuration ---

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] (%(processName)s) %(message)s',
    stream=sys.stdout,
    # also log to file 
    handlers=[
        logging.FileHandler("descjocky.log"),
        logging.StreamHandler(sys.stdout)
        ]
)

# Global, process-safe calculators
# These will be initialized once per worker process

g_mordred_calc = None
g_bluedesc_calc = None
chemopy_calc = None 

# --- Functions ---

def load_config(config_path):
    """
    Loads and validates the configuration from the .ini file.
    
    Returns:
        A dictionary of configuration settings.
    """
    if not Path(config_path).exists():
        logging.fatal(f"Configuration file not found: {config_path}")
        sys.exit(1)
    config = configparser.ConfigParser()
    config.read(config_path)
    
    settings = {}
    
    # --- [files] section ---
    try:
        files_section = config['files']
        settings['input_file'] = Path(files_section['input_file']).resolve()
        settings['mol_dir'] = Path(files_section['mol_dir']).resolve()
        
        # Optional file settings
        settings['csv_output'] = Path(files_section.get(
            'csv_output', 
            'descriptors.csv')).resolve()
        
        if not settings['input_file'].exists():
            logging.fatal(f"Input SMILES file not found: {settings['input_file']}")
            sys.exit(1)
            
    except KeyError as e:
        logging.fatal(f"Missing required configuration key: {e} in [files] section.")
        sys.exit(1)
    except Exception as e:
        logging.fatal(f"Error processing [files] section: {e}")
        sys.exit(1)

    # --- [settings] section ---
    try:
        settings_section = config['settings']
        
        # Parallelism setting
        num_workers_raw = int(settings_section.get('num_workers', 0))
        if num_workers_raw <= 0:
            settings['num_workers'] = mp.cpu_count()
        else:
            settings['num_workers'] = num_workers_raw
        
        # Do we keep intermediate files?  
        settings['remove_temp_files'] = settings_section.getboolean('remove_temp_files', False)
        
    except Exception as e:
        logging.fatal(f"Error processing [settings] section: {e}")
        sys.exit(1)

    # --- [external_tools] section ---
    settings['xtb_path'] = shutil.which('xtb')
    if settings['xtb_path'] is None:
        logging.fatal(
            "xtb executable not found in $PATH. "
        )
        sys.exit(1)

    # --- Prepare Directories ---
    settings['mol_dir'].mkdir(parents=True, exist_ok=True)
    settings['csv_output'].parent.mkdir(parents=True, exist_ok=True)
    
    logging.info(f"Configuration loaded successfully.")
    logging.info(f"Input file: {settings['input_file']}")
    logging.info(f"Output CSV: {settings['csv_output']}")
    logging.info(f"Temp dir: {settings['mol_dir']}")
    logging.info(f"Parallel workers: {settings['num_workers']}")
    return settings


def prefix_keys(d, prefix):
    """Helper function to add a prefix to all keys in a dictionary."""
    if not isinstance(d, dict):
        return {}
    return {f"{prefix}_{k}": v for k, v in d.items()}


def init_worker(): 
    """Initialize descriptor calculators for each worker process."""
    global g_mordred_calc, g_bluedesc_calc, chemopy_calc 
    g_mordred_calc = Calculator(descriptors, ignore_3D=False)
    g_bluedesc_calc = BlueDesc(ignore_3D=False)
    chemopy_calc = ChemoPy(ignore_3D=False, include_fps=True) 
    logging.info("Worker initialized with descriptor calculators.")


def run_geometry_optimization(smiles_list, config): 
    """
    Generates initial 3D conformers, writes them to an SDF file, 
    and runs a single xtb optimization for all molecules, 
    returning the path to the optimized SDF file.
    """
    
    logging.info(f"Starting geometry optimization for {len(smiles_list)} molecules.")
   
    initial_dir = config['mol_dir'] / "initial" 
    initial_dir.mkdir(parents=True, exist_ok=True)
    opt_dir = config['mol_dir'] / "optimized"
    opt_dir.mkdir(parents=True, exist_ok=True)

    successful_mols = 0 
    
    for i, smiles in enumerate(smiles_list):
        mol_id = f"mol_{i+1:04d}"
        try: 
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            if mol is None:
                logging.warning(f"Invalid SMILES: {smiles} (ID: {mol_id}). Skipping.")
                continue
            status = AllChem.EmbedMolecule(mol)
            if status != 0:
                logging.warning(f"Embedding failed for SMILES: {smiles} (ID: {mol_id}). Skipping.")
                continue
            mol.SetProp("_Name", f"{mol_id}: {smiles}")
            sdf_path = initial_dir / f"{mol_id}.sdf"
            with Chem.SDWriter(str(sdf_path)) as writer:
                writer.write(mol) 
        except Exception as e:
            logging.error(f"Error processing SMILES: {smiles} (ID: {mol_id}): {e}. Skipping.")
            # continue to next molecule 
            continue

        # Run xtb optimization 
        xtb_cmd = [ 
            config['xtb_path'], 
            str(sdf_path),
            "--opt",
            "-P", str(config['num_workers'])
        ]
        expected_opt_sdf = initial_dir / f"xtbopt.sdf"
        final_destination = opt_dir / f"{mol_id}.sdf"
        try: 
            logging.info(f"Running xtb optimization for {mol_id}.")
            result = subprocess.run(
                xtb_cmd, 
                cwd=initial_dir,
                capture_output=True,
                text=True,
                timeout=None 
            )
            if result.returncode == 0 and expected_opt_sdf.exists():
                shutil.move(str(expected_opt_sdf), str(final_destination))
                successful_mols += 1
                logging.info(f"Optimization completed for {mol_id}.")
            else:
                logging.warning(f"xtb optimization failed for {mol_id} ({smiles}): {result.stderr}")
        except Exception as e:
            logging.error(f"Error during xtb optimization for {mol_id} ({smiles}): {e}")
        finally:
            if config['remove_temp_files']:
                try:
                    shutil.rmtree(initial_dir)
                except Exception as e:
                    logging.error(f"Error removing temp files in {initial_dir}: {e}")
    if successful_mols == 0:
        logging.fatal("No molecules were successfully optimized. Exiting.")
        sys.exit(1)
    logging.info(f"Geometry optimization completed. {successful_mols} / {len(smiles_list)} molecules optimized.")
    return opt_dir

def load_optimized_molecules(opt_mol_dir_path):
    """
    Helper function to load all molecules from the optimized SDF
    directory into a list of RDKit Mol objects.
    """
    logging.info(f"Loading optimized molecules from {opt_mol_dir_path}...")
    try:
        # Glob for all SDF files and sort them numerically by their stem (e.g., '1', '2', '10')
        sdf_files = sorted(
            Path(opt_mol_dir_path).glob("*.sdf"),
            key=lambda p: int(p.stem)
        )
        
        if not sdf_files:
            logging.error(f"No .sdf files found in {opt_mol_dir_path}.")
            return []
            
        opt_mols = []
        for sdf_path in sdf_files:
            try:
                supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
                mol = supplier[0]
                if mol is not None:
                    opt_mols.append(mol)
                else:
                    logging.warning(f"Could not load molecule from {sdf_path}.")
            except Exception as e:
                logging.warning(f"Failed to process {sdf_path}: {e}")
                
        logging.info(f"Successfully loaded {len(opt_mols)} optimized molecules.")
        return opt_mols
    except Exception as e:
        logging.fatal(f"Could not read optimized molecules from {opt_mol_dir_path}: {e}")
        sys.exit(1)

def process_molecule_worker(opt_mol):
    """
    Calculates all descriptors for a single optimized molecule. 

    This function is run inside the multiprocessing pool.
    It takes one RDKit Mol object (with 3D geometry) and
    calculates all descriptors.
    """

    global g_mordred_calc, g_bluedesc_calc

    try:
        # Retrieve the original SMILES from the mol property
        smiles_string = opt_mol.GetProp("_Name")
    except (KeyError, AttributeError):
        logging.warning("Optimized mol missing '_Name' property. Re-generating SMILES.")
        try:
            smiles_string = Chem.MolToSmiles(opt_mol, isomericSmiles=True)
        except Exception:
            logging.error("Could not get SMILES from optimized mol. Skipping.")
            return None

    all_descs = {'input_smiles': smiles_string}

    try:
        # --- Calculate All Descriptors ---
        
        # PaDEL (uses SMILES, 2D/3D)
        try:
            padel_descs = padelpy.from_smiles(smiles_string, timeout=60)
            all_descs.update(prefix_keys(padel_descs, 'PaDEL'))
        except Exception as e:
            logging.warning(f"PaDEL failed for {smiles_string}: {e}")

        # Mordred (uses optimized RDKit Mol, 2D/3D)
        try:
            mordred_descs = g_mordred_calc(opt_mol).asdict()
            all_descs.update(prefix_keys(mordred_descs, 'Mordred'))
        except Exception as e:
            logging.warning(f"Mordred failed for {smiles_string}: {e}")
            
        # BlueDesc (uses optimized RDKit Mol, 2D/3D)
        try:
            # BlueDesc expects a list
            bluedesc_df = g_bluedesc_calc.calculate([opt_mol])
            if not bluedesc_df.empty:
                bluedesc_descs = bluedesc_df.iloc[0].to_dict()
                all_descs.update(prefix_keys(bluedesc_descs, 'BlueDesc'))
        except Exception as e:
            logging.warning(f"BlueDesc failed for {smiles_string}: {e}")
            
        # Pybel (uses SMILES, 2D only)
        try:
            pybel_mol = pybel.readstring('smi', smiles_string)
            pybel_descs = pybel_mol.calcdesc()
            all_descs.update(prefix_keys(pybel_descs, 'Pybel'))
        except Exception as e:
            logging.warning(f"Pybel failed for {smiles_string}: {e}")

        # Chemopy (uses optimized RDKit Mol, 2D/3D/Fingerprints)
        try:
            # expects a RDKit Mol object of molecules, outputs a Pandas DataFrame
            chemopy_descs = chemopy_calc.calculate([opt_mol])
            if not chemopy_descs.empty:
                chemopy_descs_dict = chemopy_descs.iloc[0].to_dict()
                all_descs.update(prefix_keys(chemopy_descs_dict, 'Chemopy'))
        except Exception as e:
            logging.warning(f"Chemopy failed for {smiles_string}: {e}")
            
        return all_descs

    except Exception as e:
        logging.error(f"FATAL error processing {smiles_string}: {e}", exc_info=True)
        return None

def write_csv(results, out_path):
    """
    Writes the list of descriptor dictionaries to a CSV file.
    Dynamically determines the header from all keys present.
    """
    if not results:
        logging.fatal("No results to write to CSV.")
        return 
    logging.info(f"Writing {len(results)} results to CSV: {out_path}")
    all_keys = set()
    for res in results:
        all_keys.update(res.keys())
    # Ensure a consistent order, with 'input_smiles' first
    try:
        all_keys.remove('input_smiles')
        ordered_keys = ['input_smiles'] + sorted(list(all_keys))
    except KeyError:
        ordered_keys = sorted(list(all_keys))
    try:
        with open(out_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=ordered_keys, restval='NA')
            writer.writeheader()
            writer.writerows(results)
        logging.info("CSV file written successfully.")
    except Exception as e:
        logging.fatal(f"Failed to write CSV file: {e}", exc_info=True)

def main():
    """
    Main execution function:
    1. Parse arguments
    2. Load config
    3. Run Phase 1: Geometry Optimization
    4. Run Phase 2: Descriptor Calculation
    5. Write results and cleanup
    """
    parser = argparse.ArgumentParser(description="DescJocky: A Cheminformatics Descriptor Aggregator")
    parser.add_argument(
        '-c', '--config',
        default='config.ini',
        help='Path to the configuration file (default: config.ini)'
    )
    args = parser.parse_args()
    
    # --- 1. Load config ---
    config = load_config(args.config)

    # --- 2. Read input SMILES ---
    try:
        with open(config['input_file'], 'r') as f:
            smiles_list = [line.strip() for line in f if line.strip()]
    except Exception as e:
        logging.fatal(f"Could not read input SMILES file: {e}")
        sys.exit(1)
        
    if not smiles_list:
        logging.fatal("Input SMILES file is empty.")
        sys.exit(1)
        
    logging.info(f"Found {len(smiles_list)} total molecules for processing.")

    # --- 3. Run Phase 1: Geometry Optimization ---
    try:
        opt_mol_path = run_geometry_optimization(smiles_list, config)
    except Exception as e:
        logging.fatal(f"Phase 1 (Geometry Optimization) failed: {e}", exc_info=True)
        sys.exit(1)

    # --- 4. Run Phase 2: Descriptor Calculation ---
    
    # Load all the newly optimized molecules from the single SDF
    opt_mols = load_optimized_molecules(opt_mol_path)
    
    if not opt_mols:
        logging.fatal("No molecules were successfully optimized by XTB.")
        sys.exit(1)

    logging.info(f"Starting Phase 2: Descriptor Calculation for {len(opt_mols)} molecules...")
    
    results = []
    
    try:
        # Run the descriptor calculation in parallel
        with mp.Pool(processes=config['num_workers'], initializer=init_worker) as pool:
            # map will process the list of mols in chunks and preserve order
            results = pool.map(process_molecule_worker, opt_mols)
            
    except KeyboardInterrupt:
        logging.warning("Calculation cancelled by user.")
        sys.exit(1)
    except Exception as e:
        logging.fatal(f"A critical error occurred in the processing pool: {e}", exc_info=True)
        sys.exit(1)

    # --- 5. Filter, write output, and cleanup ---
    valid_results = [r for r in results if r is not None and len(r) > 1]
    
    if not valid_results:
        logging.error("No molecules were processed successfully.")
    else:
        logging.info(
            f"Successfully processed {len(valid_results)} / {len(opt_mols)} molecules."
        )
        write_csv(valid_results, config['csv_output'])
        
    if config['remove_temp_files']:
        logging.info("Cleaning up main temporary directory...")
        try:
            shutil.rmtree(config['mol_dir'])
        except Exception as e:
            logging.error(f"Could not remove main temp dir: {config['mol_dir']}. Error: {e}")
            
    logging.info("DescJocky run complete.")


if __name__ == "__main__":
    try: 
        mp.set_start_method('spawn')
    except RuntimeError:
        pass
    print("DescJocky: A Cheminformatics Descriptor Aggregator")
    main()
