import argparse
import logging
import warnings
warnings.filterwarnings('ignore')
from helpers import *
from regions import (load_data, region_detect_mist,
                      auto_name_assign, save_ReSort,
                        plot_regions)
from stsc_utils import (run_stsc_region_level, 
                        run_stsc_finer_level)
from finer import resort_integrate

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True,
                        help="Path to count csv file or Visium output directory, with spatial/ sub-directory.")       # positional argument
    parser.add_argument('--param_fn', type=str, required=True, 
                        help='path to yaml file containing parameters for MIST and RESORT.')
    parser.add_argument('--finer', type=bool, required=False, default=False, 
                        help='whether to perform a finer cell type deconvolution or not, if yes, reference data is required.')
    parser.add_argument('--tool', type=str, required=False, default='stsc', 
                        help='base reference-based deconvolution tools to be used for RESORT.')
    parser.add_argument('--ref_count_fp', type=str, required=False, 
                        help='path to tsv file with count matrix for reference data, rows are cells, columns are genes. Required if finer == True')
    parser.add_argument('--ref_meta_fp', type=str, required=False, 
                        help='path to meta file for reference , rows: cells, column bio_celltype: cell type annotations. Required if finer == True')
    parser.add_argument('-o', '--outdir', type=str, required=True,
                        help='folder to save ReSort results')
    args = parser.parse_args()
    return args

def main():    
    args = get_args()
    setup_outdir(args.outdir, args.finer)
    rd = load_data(args.input)
    if rd is None:
        logging.info("Data paths need to be provided.")
        return
    params = fetch_params_yaml(args.param_fn)
    rd = region_detect_mist(rd, params)
    rd = auto_name_assign(rd, params['RESORT'])
    save_ReSort(rd, f"{args.outdir}/region")
    plot_regions(rd, outdir = f"{args.outdir}/region")
    
    if args.tool == 'stsc':
        run_stsc_region_level(outdir=f"{args.outdir}/region")
        if args.finer:
            run_stsc_finer_level(
                ref_count_fp = args.ref_count_fp,
                ref_meta_fp = args.ref_meta_fp,
                st_fp = f"{args.outdir}/region/ReSort_mixture_raw.tsv",
                outdir = f"{args.outdir}/finer")
            
            resort_integrate(args.outdir, params['CellHierarchy'])

if __name__ == "__main__":
    main()