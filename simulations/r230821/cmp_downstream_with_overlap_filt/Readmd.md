# Changes to pipeline

## Why
The estimated Ne trajectory using inferred IBD tend to have oscillation in the first 20
generation

## Strategy

Change the pipeline and script to filter inferred IBD segments 
that overlapped with true IBD segment with Tmrca less than 1.5.

Involved files:
- [x] env/additional_setup.md: add tskibd-filter installation method
- [x] bin/proc_dist_ne.py: need to accept true IBD files as well
- [x] main.nf: need pipe true IBD files into non-tskibd IBD callers.

## Tool: 
I have added a tool `tskibd-filter` within the `tskibd` repo:
https://github.com/bguo068/tskibd/tree/main/tskibd-filter


# run pipeline
bmibdcaller_simulations commit: 747cf6cef1ce3c84bcba3693a2926088906bab98

env: `/local/devel/bing.guo/mambaforge/envs/bmibdcaller_simulations`
run-dir: `/local/scratch/bing/bmibdcaller_simulations/simulations/r230821/cmp_downstream`
command: `nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf -profile sge -resume`

# analysis
env: `/data/bing/mambaforge/envs/py3`
run-dir: `/local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r230821/cmp_downstream_with_overlap_filt`
command: `python3 ./01_collect_info.py && python3 ./02_analyze.py`
