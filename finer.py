"""Functions to adjust external-reference-based deconvolution results using RESORT"""
import pandas as pd
import logging
from stsc_utils import load_stsc_result
from pathlib import Path

def row_norm_df(df: pd.DataFrame) -> pd.DataFrame:
    """Row normalize such that all columns sum to one for each row"""
    return df.div(df.sum(axis=1), axis=0)

def adjust_external(df_mist: pd.DataFrame, df_external: pd.DataFrame, cell_structure: dict) -> pd.DataFrame:
    """Adjust external-refernece results using RESORT's 2-step strategy."""
    cols_to_use = [c for c in df_mist.columns.tolist() if c in cell_structure]
    df_resort = []
    for c in cols_to_use:    
        sub_cts = cell_structure[c]
        df_external[list(set(cell_structure[c]) - set(df_external.columns))] = 0        
        sub_df = row_norm_df(df_external[sub_cts])
        sub_df = sub_df.apply(lambda x: x*df_mist[c], axis=0)
        df_resort.append(sub_df)
    return pd.concat(df_resort, axis=1)

def resort_integrate(out_dir: str, cell_structure: dict) -> None:
    """RESORT 2-step strategy to integrate region-level and finer cell types' results"""
    df_mist = load_stsc_result(f"{out_dir}/region")
    df_external = load_stsc_result(f"{out_dir}/finer") 
    df_integrated = adjust_external(df_mist, df_external, cell_structure)
    df_integrated.to_csv(f"{out_dir}/final/RESORT_final.tsv", sep="\t")
    logging.info(f"RESORT adjusted external-ref-based results and saved results to {out_dir}/final/")