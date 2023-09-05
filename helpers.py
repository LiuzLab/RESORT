"""Helper functions for reading parameters and set up output directories."""
import yaml
import pandas as pd
import logging
from pathlib import Path
from os.path import exists

def fetch_params_yaml(file_path: str) -> dict:
    """Read parameters from yaml file and return a dict object."""
    return yaml.safe_load(open(file_path, "r"))

def setup_outdir(outdir: str, finer=False) -> None:
    """Set up output directory by creating any parental directories."""
    Path(outdir).mkdir(parents=True, exist_ok=False)
    logging.basicConfig(filename=f'{outdir}/run.log',
                        level=logging.INFO)
    Path(f"{outdir}/region").mkdir(parents=True, exist_ok=False)
    logging.info(f'Region-level reference, results and artifacts saved at {outdir}/region_level.')
    if finer:
        Path(f"{outdir}/finer").mkdir(parents=True, exist_ok=False)
        Path(f"{outdir}/final").mkdir(parents=True, exist_ok=False)
        logging.info(f'Finer cell type deconvolution results and artifacts saved at {outdir}/region_level.')

