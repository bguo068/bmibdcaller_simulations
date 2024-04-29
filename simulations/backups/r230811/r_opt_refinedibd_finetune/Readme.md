# Goal

to further explore how LOD affect false negative rate

# json file

use json file to provide the parameter list instead of hard-coded in nextflow file

# Command

scratch_dir: /local/scratch/bing/bmibdcaller_simulations/simulations/r230811/r_opt_refinedibd_finetune

~/.local/bin/nextflow \
    /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf \
    -resume \
    -profile sge \
    -entry OPTIMIZE_REFINEDIBD \
    --refinedibd_optimize_params_json params_list.json

