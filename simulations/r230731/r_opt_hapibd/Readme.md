# bmibdcaller_simulation commit: 398eeed7b6a6b5403eaa9268abd2ad087d2b50bc
mkdir -p /local/scratch/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd
cd /local/scratch/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd
nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf  -entry OPTIMIZE_HAPIBD -profile sge -resume
rsync -aL res /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd/
