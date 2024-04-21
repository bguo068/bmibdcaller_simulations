# goal

update pipeline for benchmark IBD caller at IBD level
- allow running ishare/ibdutil within the pipeline
    - updated ishare/ibdutils so that the pop-level total ibd is also calculated 
- homogenize the optimization of all ibd callers in the `PARAM_OPTIMIZATION` subworflow
- move all caller-specific parameters to json files so that optimization can be more
flexible to explore parameter space.