#! /usr/bin/env bash
set -eEx -o pipefail

# modify this to a path that you have write permission
SCRATCHDIR=/local/scratch/bing/
PIPELINEDIR=/local/chib/toconnor_grp/bing/bmibdcaller_simulations/

source ~/conda_devel.sh
conda activate bmibdcaller_simulations

mkdir -p ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r231226  
cd ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r231226  

# generate json files
python ${PIPELINEDIR}/simulations/r231226/gen_json.py 

# workaround with https://github.com/It4innovations/hyperqueue/issues/789
# export NXF_VER=24.10.2
export NXF_VER=25.02.1-edge
# also make sure using hq version 0.17.0

## lauch hq worker via
## in a login/compute node: 
# hq server start
## in another shell of the login node
# source ~/conda_devel.sh
# conda activate bmibdcaller_simulations
# sbatch -n 30 --mem=100g --export=ALL --time 24:00:00 --wrap="hq worker start --cpus 30 --resource 'mem=sum(102400)' " # run multiple time to get more workers
# sbatch -n 128 --mem=1024000m --export=ALL --time 96:00:00 --wrap="hq worker start --cpus 128 --resource 'mem=sum(1024000)' "


# run the pipeline
nextflow  \
    ${PIPELINEDIR}/main.nf \
    -profile hq \
    -resume \
    -entry PARAM_OPTIMIZATION \
    --use_ulimit true \
    --csp_json '*.json' 

