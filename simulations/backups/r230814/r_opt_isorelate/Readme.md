
# Info for this run

bmibdcaller_simulation commit version   2855405df55156453af971d015dec7d292b59e4d
server                                  thanos (sge)
conda environment                       /local/devel/bing.guo/mambaforge/envs/bmibdcaller_simulations
nextflow command used                   nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf -resume -profile sge -entry OPTIMIZE_ISORELATE 

# Info for plotting
plot conda environment                 /data/bing/mambaforge/envs/py3
```
for imodel in 0 1 3; do 
    python3 ./compare_ibd_and_plot.py $imodel $minmac; 
    echo $imodel $minmac; 
done
```
