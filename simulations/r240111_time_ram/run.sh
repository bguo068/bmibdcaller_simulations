#! /usr/bin/env bash
set -eEx -o pipefail

# modify this to a path that you have write permission
SCRATCHDIR=/local/scratch/bing/
PIPELINEDIR=/local/chib/toconnor_grp/bing/bmibdcaller_simulations/

conda activate bmibdcaller_simulations

python3 gen_json.py

mkdir -p ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r240111_time_ram/
cd ${SCRATCHDIR}/bmibdcaller_simulations/simulations/r240111_time_ram/

# small sample size, single thread
nextflow \
    ${PIPELINEDIR}/main.nf \
    -profile clean \
    -qs 1 \
    -resume  \
    -entry WF_SP_COMPUTATION_BENCH \
    --sp_sets_json "demog_small_size.json" \
    --isorelate_minmaf 0.01 --minmaf 0.01 \
    --nchroms 3 \
    --compute_bench_large_size 0 \
    -process.cpus 1 \
    --resdir res_small_1t 

# small sample size, 10 threads
nextflow \
    ${PIPELINEDIR}/main.nf \
    -profile clean \
    -qs 1 \
    -resume  \
    -entry WF_SP_COMPUTATION_BENCH \
    --sp_sets_json "demog_small_size2.json" \
    --isorelate_minmaf 0.01 --minmaf 0.01 \
    --nchroms 3 \
    --compute_bench_large_size 0 \
    -process.cpus 10 \
    --resdir res_small_10t 

# # large sample size, 60 threads
# nextflow \
#     ${PIPELINEDIR}/main.nf \
#     -profile clean \
#     -resume  \
#     -entry WF_SP_COMPUTATION_BENCH \
#     --sp_sets_json "demog_large_size.json" \
#     --isorelate_minmaf 0.01 --minmaf 0.01 \
#     --nchroms 3 \
#     -process.cpus 64 \
#     --compute_bench_large_size 1 \
#     --resdir res_large_64t

