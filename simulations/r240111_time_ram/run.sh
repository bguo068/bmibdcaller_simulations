#! /usr/bin/env bash

set -eEx -o pipefail

source ~/conda_devel.sh
conda activate bmibdcaller_simulations

mkdir -p /local/scratch/bing/bmibdcaller_simulations/simulations/r240111_time_ram2/
cd /local/scratch/bing/bmibdcaller_simulations/simulations/r240111_time_ram2/

# NOTE: maf and isorelate_minmaf are set to 0.01 to make sure all caller use the same number of sites

# small sample size, single thread
~/.local/bin/nextflow \
    /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
    -profile clean \
    -resume  \
    -entry WF_SP_COMPUTATION_BENCH \
    --sp_sets_json "/local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240111_time_ram/demog_small_size.json" \
    --isorelate_minmaf 0.01 --minmaf 0.01 \
    --nchroms 3 \
    --compute_bench_large_size 0 \
    -process.cpus 1 \
    --resdir res_small_1t 

# small sample size, 10 threads
~/.local/bin/nextflow \
    /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
    -profile clean \
    -resume  \
    -entry WF_SP_COMPUTATION_BENCH \
    --sp_sets_json "/local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240111_time_ram/demog_small_size2.json" \
    --isorelate_minmaf 0.01 --minmaf 0.01 \
    --nchroms 3 \
    --compute_bench_large_size 0 \
    -process.cpus 10 \
    --resdir res_small_10t 

# large sample size, 60 threads
~/.local/bin/nextflow \
    /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
    -profile clean \
    -resume  \
    -entry WF_SP_COMPUTATION_BENCH \
    --sp_sets_json "/local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240111_time_ram/demog_large_size.json" \
    --isorelate_minmaf 0.01 --minmaf 0.01 \
    --nchroms 3 \
    -process.cpus 64 \
    --compute_bench_large_size 1 \
    --resdir res_large_64t


# test
# ~/.local/bin/nextflow \
#     /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
#     -profile clean \
#     -stub \
#     -resume  \
#     -entry WF_SP_COMPUTATION_BENCH \
#     --sp_sets_json "/local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240111_time_ram/demog_large_size.json" \
#     -process.cpus 1 \
#     --compute_bench_large_size 0 \
#     --resdir res_large
