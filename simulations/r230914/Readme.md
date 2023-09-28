# versions
pipeline version: 7544cc8824dfb18e981fa69df5fbbf7d76319f34
ibdutils (python) version: 

# pipeline

NOTE: pipeline is slimified to main2

to run it copy main2.nf to the root folder of the pipeline 
```
cp main2.nf ../../
```
and then 
```
cd /autofs/scratch/bing/bmibdcaller_simulations/simulations/r230914

# first
taskset -c 15-64 ~/.local/bin/nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf -resume
# results folder: opti/res

# second
taskset -c 15-64 ~/.local/bin/nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
    -resume \
    --hapibd_minseed 2.0 \
    --hapibd_minextend 1.0 \
    --hapibd_maxgap 1000 \
    --hapibd_minmarkers 100 \
# results folder: noopti/res
```

# results: see 
repo : bmibdcaller
commit: '247251729ece9fe4ae308b5ff8df1c9f21967f37'
script: astmh23/fig1.py (see plot in interactive window

