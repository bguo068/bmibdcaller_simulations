
bmibdcaler_simulation commit: dfbb43dd04d568d14c6acc2d028bf7127566495e 
```sh
mkdir -p /local/scratch/bing/bmibdcaller_simulations/simulations/r230807/r_opt_tpbwt_finetune
cd /local/scratch/bing/bmibdcaller_simulations/simulations/r230807/r_opt_tpbwt_finetune 
. ~/conda_devel.sh 
conda activate bmibdcaller_simulations
nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf -resume -profile sge -entry OPTIMIZE_TPBWT
rsync -aLP res /autofs/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r230807/r_opt_tpbwt_finetune/
cd /autofs/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r230807/r_opt_tpbwt_finetune/
conda deactivate
conda deactivate
. ~/conda_lambda.sh
conda activate py3
for i in {0..3}; do python compare_ibd_and_plot.py $i; done
```
