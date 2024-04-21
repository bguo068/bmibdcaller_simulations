#! /usr/bin/env bash

#rerun on 1/16/24 after fix mac/maf issue

. ~/conda_devel.sh
conda activate bmibdcaller_simulations
mkdir -p /local/scratch/bing/bmibdcaller_simulations/simulations/r240108/
cd /local/scratch/bing/bmibdcaller_simulations/simulations/r240108/

# Note the 2nd run of the pipeline will use the cached, shared results from first run

# optimized
~/.local/bin/nextflow \
    /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
    -resume -profile hq \
    --sp_sets_json /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240108/sp_sets.json \
    --mp_sets_json /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240108/mp_sets.json \
    --ibdne_mincm 2 \
    --filt_ibd_by_ov false \
    --resdir res_ibdnemincm2 && \
~/.local/bin/nextflow \
    /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
    -resume -profile hq \
    --sp_sets_json /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240108/sp_sets.json \
    --mp_sets_json /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240108/mp_sets.json \
    --ibdne_mincm 4 \
    --filt_ibd_by_ov false \
    --resdir res_ibdnemincm4 && \
~/.local/bin/nextflow \
    /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
    -resume -profile hq \
    --sp_sets_json /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240108/sp_sets.json \
    --mp_sets_json /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240108/mp_sets.json \
    --ibdne_mincm 2 \
    --filt_ibd_by_ov true \
    --resdir res_ibdnemincm2_filtov && \
~/.local/bin/nextflow \
    /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
    -resume -profile hq \
    --sp_sets_json /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240108/sp_sets.json \
    --mp_sets_json /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240108/mp_sets.json \
    --ibdne_mincm 4 \
    --filt_ibd_by_ov true \
    --resdir res_ibdnemincm4_filtov

# unoptimized
~/.local/bin/nextflow \
    /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
    -resume -profile hq \
    --sp_sets_json /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240108/sp_sets.json \
    --mp_sets_json /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r240108/mp_sets.json \
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


