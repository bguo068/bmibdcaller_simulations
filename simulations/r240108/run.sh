#! /usr/bin/env bash
set -eEx -o pipefail

# modify this to a path that you have write permission
SCRATCHDIR=/local/scratch/bing/
PIPELINEDIR=/local/chib/toconnor_grp/bing/bmibdcaller_simulations/

source ~/conda_devel.sh
conda activate bmibdcaller_simulations

mkdir -p ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r240108/
cd  ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r240108/

# generate json files
python3 ${PIPELINEDIR}/simulations/r240108/gen_sets_json.py

# workaround with https://github.com/It4innovations/hyperqueue/issues/789
export NXF_VER=24.10.2
# also make sure using hq version 0.17.0

## lauch hq worker via
## in a login/compute node: 
# hq server start
## in another shell of the login node
# source ~/conda_devel.sh
# conda activate bmibdcaller_simulations
# sbatch -n 30 --mem=100g --export=ALL --time 24:00:00 --wrap="hq worker start --cpus 30 --resource 'mem=sum(102400)' " # run multiple time to get more workers


# optimized 
# with different values for the mincm parameter for IBDne and 
# with or without removing segments with TMRCA < 1.5 generations
# Note the 2nd run of the pipeline will use the cached, shared results from first run
nextflow \
    ${PIPELINEDIR}/main.nf \
    -resume -profile hq \
    --sp_sets_json sp_sets.json \
    --mp_sets_json mp_sets.json \
    --ibdne_mincm 2 \
    --filt_ibd_by_ov false \
    --resdir res_ibdnemincm2 && \
nextflow \
    ${PIPELINEDIR}/main.nf \
    -resume -profile hq \
    --sp_sets_json sp_sets.json \
    --mp_sets_json mp_sets.json \
    --ibdne_mincm 4 \
    --filt_ibd_by_ov false \
    --resdir res_ibdnemincm4 && \
nextflow \
    ${PIPELINEDIR}/main.nf \
    -resume -profile hq \
    --sp_sets_json sp_sets.json \
    --mp_sets_json mp_sets.json \
    --ibdne_mincm 2 \
    --filt_ibd_by_ov true \
    --resdir res_ibdnemincm2_filtov && \
nextflow \
    ${PIPELINEDIR}/main.nf \
    -resume -profile hq \
    --sp_sets_json sp_sets.json \
    --mp_sets_json mp_sets.json \
    --ibdne_mincm 4 \
    --filt_ibd_by_ov true \
    --resdir res_ibdnemincm4_filtov

# unoptimized
nextflow \
    ${PIPELINEDIR}/main.nf \
    -resume -profile hq \
    --sp_sets_json sp_sets.json \
    --mp_sets_json mp_sets.json \
    --ibdne_mincm 2 \
    --filt_ibd_by_ov true \
    --resdir res_ibdnemincm2_filtov_unopti \
    --hapibd_minoutput 2.0 \
    --hapibd_minseed 2.0 \
    --hapibd_minextend 1.0 \
    --hapibd_maxgap 1000 \
    --hapibd_minmarkers 100 \
    --hmmibd_n 999999 \
    --hmmibd_m 5 \
    --isorelate_imiss  0.1 \
    --isorelate_vmiss  0.1 \
    --isorelate_min_snp  20 \
    --isorelate_minmaf  0.01 \
    --refinedibd_length 2.0 \
    --refinedibd_lod  3.0 \
    --refinedibd_scale  0 \
    --refinedibd_window  40.0 \
    --tpbwt_template_opts 0 \
    --tpbwt_Lm 300 \
    --tpbwt_Lf 2.0 
