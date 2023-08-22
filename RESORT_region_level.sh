source ~/envs/env_MIST/bin/activate
python3 pseudo_internal_ref.py --input $1 --param_fn $2 --outdir $3
deactivate
source ~/envs/env_stsc/bin/activat
resort_cnt=$3/ReSort_reference_raw_count.tsv
resort_label=$3/ReSort_reference_meta.tsv
st_cnt=$data_dir/ReSort_mixture_raw.tsv

stereoscope run --sc_cnt $resort_cnt --sc_labels $resort_label --st_cnt $st_cnt -o $3 -sce 5000 -n 5000 -ste 5000 -stb 100 -scb 100

