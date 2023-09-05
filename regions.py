"""Functions to load spatial transcriptomics data, and perform region-level deconvolution."""
import sys
import logging
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
sys.path.append("MIST")
from ReST import ReST

def load_data(input_path: str) -> ReST:
    """Load spatial transcriptomics data at input_path as a ReST object."""
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

def load_csv(fn: str) -> ReST:
    """Load csv-format Spatial Transcriptomics datasets, 
    index has to be x_coordinate-y_coordinate."""

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

def load_visium(indir: str) -> ReST:
    """Load Visium spatial transcriptomics datasets from SpaceRanger/out directory.
        Direcotry must contains 'spatial/' subdirectory and h5-format filtered expression
        matrix.
    """
    rd = ReST(path=indir)
    return rd

def region_detect_mist(rd: ReST, params: dict) -> ReST:
    """Process spatial transcriptomics data and detect spatial regions."""   
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

def calculate_avg_lfc(markers: list, rd: ReST, region_ind: str) -> float:
    """Calculate average markers' log fold-changes comparing 'region_ind' to all other regions."""
    marker_results = rd.region_deg_results.copy()
    marker_results = marker_results.loc[marker_results.gene.isin(markers)]
    lfc_mean = marker_results.loc[marker_results.region_ind == region_ind, 'lfc'].mean()
    return lfc_mean

def region_to_name(marker_list: dict, region: str, rd: ReST) -> str:
    """Find the best-match cell type whose marker genes enriched the most in the region"""
    best_match, max_lfc = '', 0
    for key in marker_list.keys():
        markers = marker_list[key]
        lfc = calculate_avg_lfc(markers, rd, region)
        if lfc > max_lfc:
            best_match = key
            max_lfc = lfc
    # QC: should be at least 120% enrichment
    if max_lfc < 0.26:
        best_match = 'isolated'
    return best_match

def auto_name_assign(rd: ReST, marker_list: dict):
    """Find the best-match cell type whose marker genes enriched the most in the region
        for all regions.
    """
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

def plot_regions(rd: ReST, outdir: str) -> None:
    """Plot and save detected spatial regions to output directory."""
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
    f.savefig(f'{outdir}/spatial_regions.pdf', dpi=100, bbox_inches='tight')
    plt.close()
    logging.info(f"Spatial regions plot saved at {outdir}/spatial_regions.pdf'.")

def save_ReSort(rd, outdir):
    """Save ReSort reference to output direcotry in both .csv and .tsv formats."""
    rd.save_ReSort(outdir, fmt='csv')
    rd.save_ReSort(outdir, fmt='tsv')
    logging.info("Pseudo-internal reference generated and saved.")