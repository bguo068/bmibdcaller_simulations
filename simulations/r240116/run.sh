
. ~/conda_devel.sh
conda activate bmibdcaller_simulations
mkdir -p /local/scratch/bing/bmibdcaller_simulations/simulations/r240116/
cd /local/scratch/bing/bmibdcaller_simulations/simulations/r240116/

# NOTE: to accelerate the process
# 1. use profile hq_recomb where isorelate/hmmibd were assigned more resources
# 2. use the hmmibd-rs version to allow multithreading
~/.local/bin/nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf  \
    -resume \
    -entry WF_VARY_RECOM_RATE \
    -profile hq_recomb \
    --hmmibd_version hmmibd2

# test with -stub option
# ~/.local/bin/nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf  \
#     -entry WF_VARY_RECOM_RATE \
#     -stub \
#     --hmmibd_version hmmibd2