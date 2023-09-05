"""Functions to (facilitate to) perform stereoscope-based deconvolutions."""
import logging
import subprocess
from os.path import exists
import pandas as pd
from pathlib import Path
import sys

def run_stsc(ref_count_fn: str, ref_meta_fn: str, mixture_fn: str, outdir: str) -> None:
    """Run stereoscope using subprocess."""

    if not check_valid_stsc_reference(ref_count_fn, ref_meta_fn):
        logging.info("Error: check reference file formats, indices and so on.")
        sys.exit("Error: check reference file formats, indices and so on.")

    cmd = ["stereoscope", "run", 
           "--sc_cnt", ref_count_fn,
           "--sc_labels", ref_meta_fn,
           "--st_cnt", mixture_fn,
           "--o", outdir,
           "-sce", "5000",
           "-n", "500",
           "-ste", "5000",
           "-stb", "100",
           "-scb", "100"]
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=False)
    
def run_stsc_region_level(outdir: str) -> None:
    """Run stereoscope for region-level deconvolution."""
    logging.info('Start running stsc for region-level deconvolution.')
    run_stsc(f"{outdir}/ReSort_reference_raw_count.tsv",
                f"{outdir}/ReSort_reference_meta.tsv",
                f"{outdir}/ReSort_mixture_raw.tsv", outdir
            )
    logging.info('Completed running stsc for region-level deconvolution.')

def run_stsc_finer_level(ref_count_fp: str, ref_meta_fp: str, st_fp: str, outdir: str) -> None:
    """Run stereoscope for finer-level deconvolution purely based on external reference."""
    logging.info('Start running stsc for finer-level deconvolution.')
    run_stsc(ref_count_fp,ref_meta_fp,st_fp, outdir)
    logging.info('Completed running stsc for finer-level deconvolution.')

def check_valid_stsc_reference(ref_count_path: str, ref_meta_path: str) -> bool:
    """Check if reference count path and meta path are valid."""
    def _check_tsv(fp: str) -> bool:
        """Herlper function to check if file is tab-seperated."""
        return fp.endswith("tsv") or fp.endswith("txt")
    # check format is tab-seperated
    if not (_check_tsv(ref_count_path) and _check_tsv(ref_meta_path)):
        return False
    # check file paths exist
    if not (exists(ref_count_path) and exists(ref_meta_path)):
        return False
    # check column bio_celltype exist
    meta_df = pd.read_csv(ref_meta_path, sep='\t')
    if "bio_celltype" not in meta_df.columns.tolist():
        return False
    return True
    # # check indices matching
    # count_df = pd.read_csv(ref_count_path, usecols=[0,1], sep='\t')
    # print(count_df.iloc[:3,:], meta_df.iloc[:3,:])
    # return all(count_df.iloc[:, 0] == meta_df.iloc[:, 0])

def load_stsc_result(indir: str) -> pd.DataFrame:
    """Read in stereoscope results from a directory."""
    fp = list(Path(f"{indir}/ReSort_mixture_raw/").glob("W.*.tsv"))[0]
    return pd.read_csv(fp, sep="\t", index_col=0)