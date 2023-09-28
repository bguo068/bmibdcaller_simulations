# Changes to pipeline
- no changes to the pipleline
- just add arguments to nextflow command line to temporary revert optimized parameters
to their defaults


# run pipeline
bmibdcaller_simulations commit: 7544cc8824dfb18e981fa69df5fbbf7d76319f34

env: `/local/devel/bing.guo/mambaforge/envs/bmibdcaller_simulations`
run-dir: `/local/scratch/bing/bmibdcaller_simulations/simulations/r230912/cmp_dwnstrm_ovfilt_nooptimize`
command: 
```sh
nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf \
    -profile sge -resume \
    --filt_ibd_by_ov true \
    --tpbwt_template_opts 0 \
    --tpbwt_Lm 300 \
    --tpbwt_use_phase_correction 1 \
    --hapibd_minseed 2.0 \
    --hapibd_minextend 1.0 \
    --hapibd_maxgap 1000 \
    --hapibd_minmarkers 100 \
    --refinedibd_lod 3.0 \
    --refinedibd_scale 0 \
    --refinedibd_window 40.0 \
    --hmmibd_n 999999999 \
    --hmmibd_m 5 \
    --isorelate_min_snp 20 \
    --isorelate_min_mac 20
```


# analysis
env: `/data/bing/mambaforge/envs/py3`
run-dir: `/local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r230912/cmp_dwnstrm_ovfilt_nooptimize`
command: `python3 ./01_collect_info.py && python3 ./02_analyze.py`

