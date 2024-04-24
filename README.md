# `bmibdcaller_simulations`

This repository hosts an Identity-By-Descent (IBD) caller benchmarking pipeline
written in
[Nextflow](https://github.com/nextflow-io/nextflow), 
which uses population genetic simulation and true IBD to compare the performance
of IBD segment detection methods for species with high recombination rates and
low marker densities (per genetic unit), including _Plasmodium falciparum (Pf)_.
Population genetic simulations are based on a combination of forward simulation
via
[SLiM](https://github.com/MesserLab/SLiM), 
coalescent simulation 
[msprime](https://github.com/tskit-dev/msprime) 
and their intergration via 
[pyslim](https://github.com/ tskit-dev/pyslim). 
True IBD segments are obtained from
simulated tree sequences via 
[tskibd](https://github.com/bguo068/tskibd).

The pipeline is part of a broader project `bmibdcaller`, which also includes the
following repositories:
- [bmibdcaller\_empirical](https://github.com/bguo068/bmibdcaller_empirical) 
for empirical analysis.
- [ishare/ibdutils](https://github.com/bguo068/ishare)
for effienciently comparing two sets of inferred IBD segments.
- [tskibd](https://github.com/bguo068/tskibd) 
for obtaining true ibd segments from simulated genealogical trees. 

# What are in the pipeline:

## Varying recombination rate workflow `WF_VARY_RECOM_RATE`

This workflow allows the simulation of genomes with different recombination
rates from 1e-6 to 1e-8 /bp/generation but a fixed mutation rate. The accuracy
of detected IBD segments from different IBD callers, including 
[hap-IBD](https://github.com/browning-lab/hap-ibd), 
[hmmIBD](https://github.com/glipsnort/hmmIBD), 
[isoRelate](https://github.com/bahlolab/isoRelate), 
[Refined IBD](https://faculty.washington.edu/browning/refined-ibd.html), and 
[phased IBD](https://github.com/23andMe/phasedibd), 
are evaluated by comparing them with true IBD segments to understand how
recombination rate and low per-genetic-unit marker density affect the
performance of different IBD callers.


## Parameter optimization workflow `PARAM_OPTIMIZATION`

This workflow assesses the accuracy of detected IBD segments from a given IBD
callers using different values of a single IBD caller-specific parameter or
different combinations of values of two or more IBD caller-specific parameters.
The optimal parameter value(s) are determined if these values generate IBD
segments with low and balanced error rates. This grid search approach allows us
to benchmark different IBD callers at their best status that is tuned for
genomes with high recombination rate and low marker density. 

## Computation time/RAM usage workflow `WF_SP_COMPUTATION_BENCH`

This workflow is designed to compare the computational efficiency for different
IBD detection methods (callers) by measuring runtime and peak RAM usage.


## Benchmarking workflow (default workflow)

This workflow compares the performance of different IBD callers at both
IBD-segment level (error rates/pairwise total IBD) and higher level (estimates
positive selection signal, effective population size and population structure).
It can runs with different sets of IBD calling (with optimized or nonoptimized
parameters), processing and downstream analysis parameters.

# Software environment

1. The pipeline has been tested/run under the Linux operating system. It should also
   work on MacOS. On windows, it is recommended to use WSL (Windows
   Subsystem for Linux).
2. Install conda (miniforge). See instructions from 
[here](https://github.com/conda-forge/miniforge)
3. Create a conda environment and activate the enviroment:
```sh
conda env create -f  env/bmibdcaller_simulations.yaml
conda activate bmibdcaller_simulations
```
4. Install software that is not available via conda. See notes 
[here](./env/additional_setup.md)

# How to run the pipeline

1. For assessing the impact of high recombination rate on marker density and IBD
quality: see example in
   [simulations/r240116/run.sh](./simulations/r240116/run.sh)
2. For IBD caller parameter optimization: see example in
   [simulations/r231226/run.sh](./simulations/r231226/run.sh)
3. For comparing computation efficiency: see example in
   [simulations/r240111\_time\_ram/run.sh](simulations/r240111_time_ram/run.sh)
4. For benchmarking IBD callers at both the IBD segment level and a high level (i.e.
   downstream analysis): see example in
   [simulations/r240108/run.sh](simulations/r240108/run.sh)

## Nextflow executors:
In the pipeline configuration file `nextflow.config`, several executors have been 
used to showcase how the pipeline can be run under different computing environments,
including: 
- `local` for a single-node high performance computer(HPC)
- `sge` for grid server that is managed by the Sun Grid Engine (SGE) system
- `hq` for either single-node HPC or grid server. 
   - It supports automatic worker deployment on grid server managed by PBS or slurm.
   - For other system such as SGE, it works after manual worker deployment (which I used here)
      - For example, I can run `hq server start` on the login nodes (within tmux) 
      to launch the hyperqueue server,  then run 
      `conda activate [env]; cd [workdir]` to switch to the correct working folder 
      and conda environment, and finally launch as many worker jobs as needed by running
      `qsub -b yes -cwd -V -l h='computing_node_name' hq worker start `.
   - More document can be found 
   [here](https://it4innovations.github.io/hyperqueue/v0.18.0/deployment/worker/#deploying-a-worker-using-pbsslurm)

# Examining result folders

The path of output folders for each process can be found by checking the
`publishDir` directive under each process in the `main.nf` file. The file name
of the `output` can be found in the `path` item(s) of the output channel(s).

For instance,
the `publishDir` directive and `output` for the `RUN_IBDNE` process is: 
```
    publishDir "${resdir}/${args.genome_set_id}_${label}/${ibdcaller}/ne_output/"
    ...
    output:
        tuple val(label), val(ibdcaller), val(are_peaks_removed), path("*.ne")
```
which may correspond to the following result file
```
simulations/r240108/res_ibdnemincm2/10000_sp_neu/hmmibd/ne_output/10000_2.0_10.0_none_orig.ne

 \--------------------------------/
           resdir
                                    \---/
                                genome_set_id

                                          \----/
                                          label

                                                 \----/
                                               ibdcaller
                                                                  \-------------------------/
                                                                       output/path item
```

In particular, some of them are important for further downstream analyses or
summarization:
- Effective population size results
   - `*.ne` files under (sub)folders named `ne_output`
- InfoMap community detection results
   - `*.pq` parquet files under (sub)folder named `ifm_output`
- Selection scanning
   - can be obtained from analyzing `*.ibdobj` files under (sub)folders named
   `ifm_input` or `ibdne_ibd`, and `*.vcf.gz` files under (subfolders) named
   `vcf` using the python package
   [`ibdutils`](https://github.com/bguo068/ibdutils)
- IBD segments error rates (overlapping rates),  Pairwise total IBD, and
popluation-level total IBD per given length ranges
   - overlapping rates: `*.ovcsv` files under subfolders named `cmpibd`
   - pairwise total IBD: `*.pairtotibdpq` file under subfolder named `cmpibd`
   - length-specific population-level total IBD: `*.poptotibdcsv` file under
   subfolder named `cmpibd`
- Run time and memory usage results
   - `time_output.txt` files under subfolders of parent folder named `ibdcalltime`.
