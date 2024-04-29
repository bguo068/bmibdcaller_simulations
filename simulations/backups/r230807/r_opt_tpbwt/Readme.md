
bmibdcaler_simulation commit: 9c3cca109c2ffa02ba76fe27b379f330fe1aeff5
```sh
mkdir -p /local/scratch/bing/bmibdcaller_simulations/simulations/r230807/r_opt_tpbwt
cd /local/scratch/bing/bmibdcaller_simulations/simulations/r230807/r_opt_tpbwt 
. ~/conda_devel.sh 
conda activate bmibdcaller_simulations
nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf -resume -profile sge -entry OPTIMIZE_TPBWT
```
