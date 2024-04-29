# bmibdcaller_simulation commit: c52d6f642450b78165cd0b45db9201023ddd0434
mkdir -p /local/scratch/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd
cd /local/scratch/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd
nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf  -entry OPTIMIZE_HAPIBD -profile sge -resume
rsync -aL res /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd/
