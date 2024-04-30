#! /usr/bin/env bash
set -eEx -o pipefail

# modify this to a path that you have write permission
SCRATCHDIR=/local/scratch/bing/
PIPELINEDIR=/local/chib/toconnor_grp/bing/bmibdcaller_simulations/

conda activate bmibdcaller_simulations

mkdir -p ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r231226  
cd ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r231226  

# generate json files
python gen_json.py 

# run the pipeline
nextflow  \
    ${PIPELINEDIR}/main.nf \
    -profile hq \
    -resume \
    -entry PARAM_OPTIMIZATION \
    --csp_json '*.json' 
