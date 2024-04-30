#! /usr/bin/env bash
set -eEx -o pipefail

conda activate bmibdcaller_simulations

# modify this to a path that you have write permission
SCRATCHDIR=/local/scratch/bing/
PIPELINEDIR=/local/chib/toconnor_grp/bing/


mkdir -p ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r240116/
cd ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r240116/

# NOTE: to accelerate the process use profile hq_recomb where isorelate is
# assigned with more resources
nextflow ${PIPELINEDIR}/main.nf  \
    -resume \
    -entry WF_VARY_RECOM_RATE \
    -profile hq_recomb
