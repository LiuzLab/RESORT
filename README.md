# RESORT
An approach to enhance reference-based spatial deconvolution methods

## Reference

1. Wang, Linhua, et al. "Accurate cell type deconvolution in spatial transcriptomics using a batch effect-free strategy." bioRxiv (2022): 2022-12.

2. Wang, L., Maletic-Savatic, M. & Liu, Z. Region-specific denoising identifies spatial co-expression patterns and intra-tissue heterogeneity in spatially resolved transcriptomics data. Nat Commun 13, 6912 (2022). https://doi.org/10.1038/s41467-022-34567-0

## Dependencies
Make sure python@3.8.10 is used. 

To automatic ReSort's pipeline, we will install a virtual environment:

```console
$ RESORT_dir=$(pwd)
$ python3 -m venv env_resort
$ source env_MIST/bin/activate
(env_resort)$
```

Make sure `pip3` is using the virtual environment:
```console
(env_resort)$ which pip3
>>> env_resort/bin/pip3
```

Install the dependencies from `requirements.txt` by running `pip3 install -r $RESORT_dir/requirements.txt`.

To test if MIST dependencies are installed:

```console
(env_resort)$ cd $RESORT_dir/
(env_resort)$ python
Python 3.8.10 (default) 
[GCC 9.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import sys
>>> sys.path.append("MIST/")
>>> from ReST import ReST
>>> quit()
```

### Install Stereoscope (stsc)

To run RESORT with the latest stereoscope version, please clone `stsc` into a directory `/path/to/stsc` by running `git clone https://github.com/almaan/stereoscope.git`. To install `stereoscope`:

```console
(env_resort)$ cd /path/to/stsc
(env_resort)$ ./setup.py install
```

To test stereoscope is installed: 

```console
(env_resort)$ python -c "import stsc; print(stsc.__version__)"
stereoscope : 0.3.1
(env_resort)$ stereoscope test
successfully installed stereoscope CLI
```

For any issues occured during installation of stereoscope, please refer to [stsc-issues](https://github.com/almaan/stereoscope/issues) or open an issue in RESORT's issue panel.


### Other R programming based deconvolution methods' dependencies

Although we provided adaptors to integrate ReSort with **R**-based deconvolution methods, the dependencies for each of the methods in `deconv_models/run*.R` need to be installed separately from their sources. Please refer to the home pages for each method below to install the dependencies.

* [MuSiC] (https://xuranw.github.io/MuSiC/articles/MuSiC.html)
* [RCTD/spacexr] (https://github.com/dmcable/spacexr)
* [Spotlight] (https://marcelosua.github.io/SPOTlight/)
* [SpatialDWLS/Giotto] (https://rubd.github.io/Giotto_site/)

<!-- ## Input data format
1. For 10X Visium, Space Ranger `Folder` with the following contents:
  - [Folder]/spatial/tissue_positions_list.csv
  - [Folder]/filtered_feature_bc_matrix.h5
2. General spatial transcriptomics data:
  - counts.csv - gene expression data frame in Pandas.DataFrame format, where each row is a spot and each column is a gene. Spots' indices are in the format of `axb` where `a` is `x coordinate` and `b` is `y coordinate`. -->

## Run RESORT

To run RESORT, call:

```
  python RESORT.py --input INPUT_PATH \
                                --outdir OUTPUT_DIR \
                                --param_fn PARAM_YAML_FN \
                                --finer OPTIONAL \
                                --ref_count_fp OPTIONAL \
                                --ref_meta_fp OPTIONAL
```

Parameters:

    * input: path to the input data, either of the following two types:

      1. For 10X Visium, Space Ranger `Folder` with the following contents:
          - [Folder]/spatial/tissue_positions_list.csv
          - [Folder]/filtered_feature_bc_matrix.h5
      2. General spatial transcriptomics data:
        - counts.csv - gene expression data frame in Pandas.DataFrame format, where each row is a spot and each column is a gene. Spots' indices are in the format of `axb` where `a` is `x coordinate` and `b` is `y coordinate`.
    
    * outdir: path to direcotry for saving artifacts, results, and plots.
    * tool: reference-based deconvolution algorithm used by RESORT, default as `stsc`.
    * param_fn: path to yaml file containing parameters for `MIST``, `RESORT` region marker list and `Cell type tree`.
    * finer: True of False, whether performing finer cell type deconvolution or not, if yes, an external reference is required.
    * ref_count_fp: path to reference count matrix in tsv format, required if finer is True. rows: cells, columns: genes.
    * ref_meta_fp: path to reference meta matrix in tsv format, required if finer is True. rows: cells, column `bio_celltype`: cell type annotation.

### Region-level reference-free deconvolution

To run reference-free region-level deconvolution:

```console
  python RESORT.py --input INPUT_PATH \
                                --outdir OUTPUT_DIR \
                                --param_fn PARAM_YAML_FN
```

An yaml file is required to specify MIST and RESORT parameters, please refer to `params.yaml` as a template to tune the parameters.
This yaml file contains three fields. Field `MIST` and `RESORT` are required for region-level deconvoution.
The field `MIST` spcifies parameters related to the package [MIST]() for region detection.

```console
MIST:
  species: 'Human' #either Mouse or Human
  n_pcs: 5 # number of PCs to be used
  min_sim: 0.5 # minimum of similarity scores for thresholding
  max_sim: 0.91 # maximum of similarity scores for thresholding
  min_cell_count: 3
  min_read_count: 200
  min_region_size: 20 # minimum region size
  hvg_prop: 0.9 #top 90% of hvgs used
  step_size: 0.02 #search step of the weights
```

The field `RESORT` specifies marker genes for each of the expected region-level cell types from the users' input:

```
RESORT:
  Cancer: ['TM4SF1']
  Ductal: ['CRISP3']
  Other Normal: ['PRSS1']
```

RESORT will require these marker genes to annotate each of the detected regions. 

### Finer cell type deconvolution

To run finer cell type deconvolution, a single-cell RNA-seq data and the cell type annotation file are needed as an external reference. Please refer to [stsc]() for detailed formats of such files. In brief, 
  * count.tsv: rows are cells, columns are genes, whose symbols match with the ST data.
  * meta.tsv: rows are cells and match with the `count.tsv`, one column must be `bio_celltype` denoting cell type.

The following command will yield the final results for finer cell type deconvolution on the PDAC-A ST sample:

```console
  python RESORT.py --input data/PDACA/PDACA_ST.csv \
                    --outdir results/PDACA/ \
                    --param_fn params.yaml \
                    --finer True \
                    --ref_count_fp data/PDACA/reference/pdacb_count.tsv \
                    --ref_meta_fp data/PDACA/reference/pdacb_meta.tsv
```

In finer cell type deconvolution, another field is required in the `params.yaml`, which is `CellHierarchy`:

```
CellHierarchy: 
  Cancer: ['Cancer']
  Ductal: ['Ductal - APOL1 high', 'Ductal - CRISP3 high', 'Ductal - MHC Class II', 'Ductal - terminal ductal like']
  Other Normal: ['Acinar cells', 'Endocrine cells', 'Endothelial cells', 'Macrophages', 'Mast cells', 'Monocytes', 'T cells & NK cells', 'Tuft cells', 'mDCs', 'pDCs']
```

It specifies the hierarchy of the cell types among the region-level cell types and the finer cell types. All the finer cell types appeared in the cell hierarchy should also exist in the annotation file for the external reference. 

The results will be saved to `results/PDACA/final` in this case.