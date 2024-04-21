
source ~/conda_devel.sh 
conda activate bmibdcaller_simulations
mkdir -p /local/scratch/bing/bmibdcaller_simulations/simulations/r231226  
cd /local/scratch/bing/bmibdcaller_simulations/simulations/r231226  

# generate json files
python gen_json.py 

# run the pipeline
~/.local/bin/nextflow  /local/chib/toconnor_grp/bing/bmibdcaller_simulations/main2.nf -profile hq -resume -entry PARAM_OPTIMIZATION --csp_json '*.json' 

# NOTE: 
# check whether using rare variants improve IBD quality for hmmIBD: hmmibd_origalgmaf001.json
