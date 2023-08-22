import sys
sys.path.append("MIST")
from ReST import ReST
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import yaml
import argparse
import logging

def get_yaml(file_path):
    return yaml.safe_load(open(file_path, "r"))

def get_args(file_path):
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='', required=True)           # positional argument
    # parser.add_argument('--visium', type=str, default='', required=False)
    parser.add_argument('--param_fn', type=str, required=True, help='path to yaml file containing parameters for MIST and RESORT.')
    parser.add_argument('-o', '--outdir', type=str, required=True,
                        help='folder to save ReSort results')
    args = parser.parse_args()
    return args

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
    return rd

def plot_regions(rd):
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

def load_data(args):
    input_path = args.input
    if input_path.endswith('.csv'):
        return load_csv(args.csv)
    else:
        try: 
            return load_visium(args.visium)
        except:
            print("Wrong input data.")
            return False

def main(args):
    rd = load_data(args)
    if not rd:
        logging.info("Data paths need to be provided.")
        return
    logging.info("Data successfully loaded.")
    params = get_yaml(args.param_fn)
    rd = region_detect_mist(rd, params)
    logging.info("Regions detected.")
    rd = auto_name_assign(rd, params['RESORT'])
    logging.info("Regions' names matched.")
    save_ReSort(rd, args.outdir)
    logging.info("Pseudo-internal reference generated and saved.")
    f = plot_regions(rd)
    f.savefig(f'{args.outdir}/region-levels.pdf', dpi=100, bbox_to_inches='tight')
    logging.info(f'Reference and artifacts saved at {args.outdir}/.')

if __name__ == "__main__":
    args = get_args()
    logging.basicConfig(filename=f'{args.outdir}/pseudo_internal_ref.log',level=logging.info)
    main(args)