# Changes to pipeline
1. add/update params to pipeline
    - update pipeline-level parameter default to optimized parameters
    - updated files: `main.nf` and `bin/call_isorelate.py`
2. make sure downstream analysis steps/script/code consistent with posseleff 
manuscript
    - key updated files
```
proc_dist_ne.py
proc_infomap.py
run_infomap.py
```

# run pipeline
bmibdcaller_simulations commit: 8ffe4d17f7eb59b1193c38a8303e826c9c627e44

env: `/local/devel/bing.guo/mambaforge/envs/bmibdcaller_simulations`
run-dir: `/local/scratch/bing/bmibdcaller_simulations/simulations/r230821/cmp_downstream`
command: `nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf -profile sge -resume`

# analysis
env: `/data/bing/mambaforge/envs/py3`
run-dir: `/local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r230821/cmp_downstream_no_overlap_filt`
command: `python3 ./01_collect_info.py && python3 ./02_analyze.py`
