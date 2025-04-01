#! /usr/bin/env bash
set -eEx -o pipefail

# activate conda environemtn
# source ~/conda_devel.sh
conda activate bmibdcaller_simulations

# modify this to a path that you have write permission
SCRATCHDIR=/local/scratch/bing/
PIPELINEDIR=/local/chib/toconnor_grp/bing/bmibdcaller_simulations


mkdir -p ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r240116/
cd ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r240116/

# NOTE: to accelerate the process use profile hq_recomb where isorelate is
# assigned with more resources
# Also use hmmibd-rs (here named hmmibd2) to accerlate IBD calling for smaller recombination cases
# https://github.com/bguo068/hmmibd-rs 
export NXF_VER=25.02.1-edge
nextflow ${PIPELINEDIR}/main.nf  \
    -resume \
    -entry WF_VARY_RECOM_RATE \
    -profile hq_recomb \
    --hmmibd_version hmmibd2
