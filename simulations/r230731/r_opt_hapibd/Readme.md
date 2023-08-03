# bmibdcaller_simulation commit: e13e7e91fb684ec7cd7bdbfe334877c6c440f3f5
mkdir -p /local/scratch/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd
cd /local/scratch/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd
nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf  -entry OPTIMIZE_HAPIBD -profile sge -resume

rsync -aL res /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd/
