# RESORT
An approach to enhance reference-based spatial deconvolution methods

## Reference

1. Wang, Linhua, et al. "Accurate cell type deconvolution in spatial transcriptomics using a batch effect-free strategy." bioRxiv (2022): 2022-12.

2. Wang, L., Maletic-Savatic, M. & Liu, Z. Region-specific denoising identifies spatial co-expression patterns and intra-tissue heterogeneity in spatially resolved transcriptomics data. Nat Commun 13, 6912 (2022). https://doi.org/10.1038/s41467-022-34567-0

## Dependencies

### MIST dependencies

We recommend using an virtual environment for MIST's dependencies

```console
$ cd MIST/
$ python3 -m venv env_MIST
$ source env_MIST/bin/activate
(env_MIST)$ pip install -r requirements.txt
```

To test if MIST dependencies are installed:

```console
(env_MIST)$ python
[GCC 9.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from ReST import ReST
>>>
```

### Stereoscope dependencies

We recommend using an virtual environment for stereoscope installation. To run RESORT with the latest stereoscope version, please clone `stsc` into home directory `~/stsc` or other user customized folder by running `git clone https://github.com/almaan/stereoscope.git`. To install `stereoscope`:

```console
$ cd ~/stsc/
$ python3 -m venv (env_stsc)
(env_stsc)$ source env_stsc/bin/activate
(env_stsc)$ ./setup.py install
```

To test stereoscope is installed: 

```console
(env_stsc)$ python -c "import stsc; print(stsc.__version__)"
stereoscope : 0.3.1
(env_stsc)$ stereoscope test
successfully installed stereoscope CLI
```

For any issues occured during installation of stereoscope, please refer to [stsc-issues](https://github.com/almaan/stereoscope/issues) or open an issue in RESORT's issue panel.


### Other R programming based deconvolution methods' dependencies

Although we provided adaptors to integrate ReSort with **R**-based deconvolution methods, the dependencies for each of the methods in `deconv_models/run*.R` need to be installed separately from their sources. Please refer to the home pages for each method below to install the dependencies.

* [MuSiC] (https://xuranw.github.io/MuSiC/articles/MuSiC.html)
* [RCTD/spacexr] (https://github.com/dmcable/spacexr)
* [Spotlight] (https://marcelosua.github.io/SPOTlight/)
* [SpatialDWLS/Giotto] (https://rubd.github.io/Giotto_site/)

### Input data format
1. For 10X Visium, Space Ranger `Folder` with the following contents:
  - [Folder]/spatial/tissue_positions_list.csv
  - [Folder]/filtered_feature_bc_matrix.h5
2. General spatial transcriptomics data:
  - counts.csv - gene expression data frame in Pandas.DataFrame format, where each row is a spot and each column is a gene. Spots' indices are in the format of `axb` where `a` is `x coordinate` and `b` is `y coordinate`.
