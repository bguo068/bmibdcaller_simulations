# bmibdcaller_simulation commit:9430a1cfe214818cdbae16d9394b550b9eba45d5
mkdir -p /local/scratch/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd_finetune_minmarker
cd /local/scratch/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd_finetune_minmarker
nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf  -entry OPTIMIZE_HAPIBD -profile sge -resume
rsync -aL res /local/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r230731/r_opt_hapibd_finetune_minmarker/
