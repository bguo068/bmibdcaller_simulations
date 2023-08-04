# Questions

1. Are optimized parameters demographic model-specific? Would parameters
optimized for single pop module be different from those for multiple pop model.
    - What are the NP/FP rate if we the simulation model from hap-ibd paper?
    - url: https://www.sciencedirect.com/science/article/pii/S0002929720300525

```
1. With uk populiaton model and human evolutionary parameters, hapibd has good quality
2. With single/multiple population model, hapidbd results highly depend on minmarker parameter
3. maf > 0.1 generate worse results than maf>0.01 
    - maf >= 0.001 (a.k.a no filter) can generate worse result, for instance (data not included in git)
    in the multiple population model 
4. gene conversion versus recombination only does not change the results much (in uk human model)
5. for pf single/mutiple model, minmarker = 70 or 80 seems to the best parameter.

NOTE:
"./r_opt_hapibd" supports points 1/2/3/4
"./r_opt_hapibd_finetune_minmarker" supports point 5
```

2. What parameters should be optimized for each of the to-be-compared IBD
callers, including hmmibd, refined-ibd, isoRelate, and tpbwt.
    - for hmmibd
        - url: https://github.com/glipsnort/hmmIBD
        - how `n` change the results
        - how `m` change the results
    - refined-ibd
        - url: https://faculty.washington.edu/browning/refined-ibd/refined-ibd.04Dec18.pdf
        - `scale` parameter (for haplotype frequency model)
        - `lod`: minimum LOD score for reportd IBD segments
    - isoRelate 
        - url: https://github.com/bahlolab/isoRelate/blob/master/vignettes/introduction.Rmd
        - `minimum.snps` 
        - `minimum.length.bp`
        - `error`
    - TPBWT:
        - url: https://github.com/23andMe/phasedibd
        - `use_phase_correction`
        - `Lm`:  minimum number of SNPs a matching subsegment must span to be included in an IBD segment.
        - `Lf`: The minimum genetic length (cM) of an IBD segment.
3. What is the relationship between false positive rate and sensitivity in
selection detection?
    - positive selection simulation
    - how to measure sensitivity in selection detection?
    - relationship of false negative/false positive rates and whether is detectable
4. How does IBD segments breakdown or merging affect Ne estimates and other time-specific analyses?
    - Define and calculate the rate of IBD segment breakdown and IBD merging
    - how are these rate related to bias in Ne estimates and IBD length - time
    relationship? What is the metric to analyze this?
5. Merging IBD segments can cause overestimating estimates of relatedness.
    - isorelate infer IBD in monclonal infection vs polyclonal infection
    simulations
6. hmmibd model suffers from deconvolution errors.
    - True dominant clone and inferred dominant clone
