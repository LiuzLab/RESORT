import sys
import yaml
import pandas as pd
import seaborn as sns
import argparse
import logging
import warnings
warnings.filterwarnings('ignore')
import subprocess
from pathlib import Path
from matplotlib import pyplot as plt
from os.path import exists
from shutil import rmtree
sys.path.append("MIST")
from ReST import ReST

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

def get_yaml(file_path):
    return yaml.safe_load(open(file_path, "r"))

def load_csv(fn: str):
    df = pd.read_csv(fn, index_col=0)
    df.columns = [c.split()[0] for c in df.columns]
    xs, ys = [],[]
    for ind in df.index.tolist():
        xs.append(int(ind.split('x')[0]))
        ys.append(int(ind.split('x')[1]))
    meta_df = pd.DataFrame({'array_col': xs, "array_row": ys}, index = df.index)
    gene_df = pd.DataFrame({'gene': df.columns.tolist()}, index=df.columns)
    rd = ReST(counts=df, coordinates=meta_df, gene_df=gene_df)
    return rd

def load_visium(folder: str):
    rd = ReST(path=folder)
    return rd

def save_ReSort(rd, folder):
    rd.save_ReSort(folder, fmt='csv')
    rd.save_ReSort(folder, fmt='tsv')
    logging.info("Pseudo-internal reference generated and saved.")

def region_detect_mist(rd, params):   
    rd.preprocess(species=params['MIST']['species'],
                   hvg_prop=params['MIST']['hvg_prop'], 
                   n_pcs=params['MIST']['n_pcs'],
                    min_cell_count=params['MIST']['min_cell_count'],
                      min_read_count=params['MIST']['min_read_count'])

    rd.extract_regions(min_sim = params['MIST']['min_sim'],
                        min_size=params['MIST']['min_region_size'],
                          gap=params['MIST']['step_size'])
    
    rd.assign_region_colors()
    rd.extract_regional_markers()
    logging.info("Regions detected, markers extracted.")
    return rd

def calculate_avg_lfc(markers, rd, region_ind):
    marker_results = rd.region_deg_results.copy()
    marker_results = marker_results.loc[marker_results.gene.isin(markers)]
    lfc_mean = marker_results.loc[marker_results.region_ind == region_ind, 'lfc'].mean()
    return lfc_mean

def region_to_name(marker_list, region, rd):
    best_match, max_lfc = '', 0
    for key in marker_list.keys():
        markers = marker_list[key]
        lfc = calculate_avg_lfc(markers, rd, region)
        if lfc > max_lfc:
            best_match = key
            max_lfc = lfc
    if max_lfc < 0.26:
        best_match = 'isolated'
    return best_match

def auto_name_assign(rd, marker_list):
    region_match_name = {"isolated": "isolated"}
    region_inds = list(set(rd.adata.obs.region_ind))
    name_dict = {}

    for region_ind in region_inds:
        if region_ind not in region_match_name:
            best_match = region_to_name(marker_list, region_ind, rd)
            if best_match in name_dict.keys():
                matched_name = f'{best_match}_{name_dict[region_ind]}'
                name_dict[region_ind] += 1
            else:
                matched_name = best_match
                name_dict[region_ind] = 1
            region_match_name[region_ind] = matched_name
            print(f"Region index {region_ind} named as {matched_name}.")

    rd.manual_assign_region_names(region_match_name)
    rd.adata.obs['regions_new'] = rd.adata.obs.region_ind.map(region_match_name)
    logging.info("Regions' names matched.")
    return rd

def plot_regions(rd, outdir):
    f, axs = plt.subplots(1,2, figsize=(11,5))
    region_df = rd.adata.obs.copy()

    sns.scatterplot(data=region_df, x="array_col", y="array_row", hue='region_ind',
                    palette='tab20',
                    ax=axs[0])

    sns.scatterplot(data=region_df, x="array_col", y="array_row", hue='regions_new',
                    palette="tab10", ax=axs[1])

    axs[0].axis("off")
    axs[1].axis("off")

    axs[0].legend(frameon=True, fontsize=12, bbox_to_anchor=(1.02, 0.8))
    axs[0].set_title("MIST detected regions")

    axs[1].legend(frameon=True, fontsize=12, bbox_to_anchor=(1.02, 0.8))
    axs[1].set_title("ReSort matched regions")
    f.savefig(f'{outdir}/region-levels.pdf', dpi=100, bbox_inches='tight')
    plt.close()

def load_data(input_path: str):
    if input_path.endswith('.csv'):
        rd = load_csv(input_path)
    else:
        try: 
            rd = load_visium(input_path)
        except:
            print("Wrong input data.")
            return None
    logging.info("Data successfully loaded.")
    return rd

def setup_outdir(outdir, finer=False):
    if exists(outdir):
        rmtree(outdir)
    Path(outdir).mkdir(parents=True, exist_ok=False)
    logging.basicConfig(filename=f'{outdir}/run.log',
                        level=logging.INFO)
    Path(f"{outdir}/region_level").mkdir(parents=True, exist_ok=False)
    logging.info(f'Region-level reference, results and artifacts saved at {outdir}/region_level.')
    if finer:
        Path(f"{outdir}/finer_level").mkdir(parents=True, exist_ok=False)
        logging.info(f'Finer cell type deconvolution results and artifacts saved at {outdir}/region_level.')

def run_stsc(ref_count_fn, ref_meta_fn, mixture_fn, outdir):
    if not check_valid_reference(ref_count_fn, ref_meta_fn):
        logging.info("Error: check reference file formats, indices and so on.")
        return
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
    
def run_stsc_region_level(outdir):
    logging.info('Start running stsc for region-level deconvolution.')
    run_stsc(f"{outdir}/ReSort_reference_raw_count.tsv",
                f"{outdir}/ReSort_reference_meta.tsv",
                f"{outdir}/ReSort_mixture_raw.tsv", outdir
            )
    logging.info('Completed running stsc for region-level deconvolution.')

def run_stsc_finer_level(ref_count_fp, ref_meta_fp, st_fp, outdir):
    logging.info('Start running stsc for finer-level deconvolution.')
    run_stsc(ref_count_fp,ref_meta_fp,st_fp, outdir)
    logging.info('Completed running stsc for finer-level deconvolution.')
    return

def check_tsv(fp: str) -> bool:
    return fp.endswith("tsv") or fp.endswith("txt")

def check_valid_reference(ref_count_path: str, ref_meta_path: str):
    # check format is tab-seperated
    if not (check_tsv(ref_count_path) and check_tsv(ref_meta_path)):
        return False
    # check file paths exist
    if not (exists(ref_count_path) and exists(ref_meta_path)):
        return False
    # check column bio_celltype exist
    meta_df = pd.read_csv(ref_meta_path)
    if "bio_celltype" not in meta_df.columns.tolist():
        return False
    # check indices matching
    count_df = pd.read_csv(ref_count_path, usecols=[0,1])
    return all(count_df.iloc[:, 0] == meta_df.iloc[:, 0])

# todo
def resort_integrate():
    return

def main():    
    args = get_args()
    setup_outdir(args.outdir, args.finer)
    rd = load_data(args.input)
    if rd is None:
        logging.info("Data paths need to be provided.")
        return
    params = get_yaml(args.param_fn)
    rd = region_detect_mist(rd, params)
    rd = auto_name_assign(rd, params['RESORT'])
    save_ReSort(rd, f"{args.outdir}/region_level")
    plot_regions(rd, outdir = f"{args.outdir}/region_level")
    
    if args.tool == 'stsc':
        run_stsc_region_level(outdir=args.outdir)
        if args.finer:
            run_stsc_finer_level(
                ref_count_fp = args.ref_count_fp,
                ref_meta_fp = args.ref_meta_fp,
                st_fp = f"{args.outdir}/regioin_level/ReSort_mixture_raw.tsv",
                out_dir = f"{args.outdir}/finer_level")
            
            #resort_integrate(args.outdir, params['Cell_structure'])


if __name__ == "__main__":
    main()