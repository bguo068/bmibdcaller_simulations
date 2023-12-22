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
```

# first
taskset -c 15-64 ~/.local/bin/nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf -resume
# results folder: opti/res

# updated on Dec 21. The new run will run all IBD callers (NOT only hmmibd and hapibd)
- files changes:
	- main2.nf (to uncomment gneome sets and IBD callers
	- res dir is updated 
		- old version: `res_old`
		- new version: `res`

```
taskset -c 15-64 ~/.local/bin/nextflow /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf \
    -resume \
    --tpbwt_template_opts 0 \
    --tpbwt_Lm 300 \
    --tpbwt_use_phase_correction 1 \
    --hapibd_minseed 2.0 \
    --hapibd_minextend 1.0 \
    --hapibd_maxgap 1000 \
    --hapibd_minmarkers 100 \
    --refinedibd_lod 3.0 \
    --refinedibd_scale 0 \
    --refinedibd_window 40.0 \
    --hmmibd_n 999999999 \
    --hmmibd_m 5 \
    --isorelate_min_snp 20 \
    --isorelate_min_mac 20
```
# results folder: noopti/res

# results: see 

repo : bmibdcaller
commit: '247251729ece9fe4ae308b5ff8df1c9f21967f37'
script: astmh23/fig1.py (see plot in interactive window)


# results: 

Note: the genome set with the smallest recombination rate is not run this time
as the genome size in bp is too large for run slow IBD caller, especially
isorelate.

Note2: even with the lowest recombination rate group comment out. Isorelate
still is too slow for the second smallest recombination rate group. See below
summary for how many chromosomes are done for each set.

```
Gsid   Caller        Count
40002  hapibd        14
       hmmibd        14
       isorelate      8
       refinedibd    14
       tpbwtibd      14
       tskibd        14
40003  hapibd        14
       hmmibd        14
       isorelate     12
       refinedibd    14
       tpbwtibd      14
       tskibd        14
40004  hapibd        14
       hmmibd        14
       isorelate     14
       refinedibd    14
       tpbwtibd      14
       tskibd        14
40005  hapibd        14
       hmmibd        14
       isorelate     14
       refinedibd    14
       tpbwtibd      14
       tskibd        14
40006  hapibd        14
       hmmibd        14
       isorelate     14
       refinedibd    14
       tpbwtibd      14
       tskibd        14
40007  hapibd        14
       hmmibd        14
       isorelate     14
       refinedibd    14
       tpbwtibd      14
       tskibd        14
```
