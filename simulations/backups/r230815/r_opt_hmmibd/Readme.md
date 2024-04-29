
# Info for this run

bmibdcaller_simulation commit version   adff16d8581d1077fc817096f6583b9f68215329
server                                  lambda
conda environment                       /local/devel/bing.guo/mambaforge/envs/bmibdcaller_simulations
nextflow command used                   nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main.nf -resume -entry OPTIMIZE_HMMIBD 

# Info for plotting
```
for imodel in 0 1 3; do 
    for minmac in 1 20 200; do 
        ./compare_ibd_and_plot.py $imodel $minmac; 
        echo $imodel $minmac; 
    done; 
done
```
