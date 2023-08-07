#! /usr/bin/env python3

import phasedibd
import pandas as pd
import numpy as np
import allel
import argparse
from pathlib import Path
import sys
import subprocess
from pathlib import Path

# --------------------- parse arguments -----------------------------
parser = argparse.ArgumentParser(
    """
call_tpbwt.py does to three things
1. automatically generate the genetic map
2. call tpbwt on a single-chromosome vcf file
3. reformat ibd file to a similar format as used by tree-based raw ibd, fake columns are
    - df_ibd["Ancestor"] = 99999
    - df_ibd["Tmrca"] = 100
    - df_ibd["HasMutation"] = 0
The output will be ${chrno}.ibd, ${chrno}.map and ${chrno}.log file
"""
)
parser.add_argument("--vcf", type=str, required=True, help="input vcf file")
parser.add_argument("--r", type=float, required=True)
parser.add_argument("--seqlen", type=int, required=True)
parser.add_argument("--chrno", type=int, required=True)
parser.add_argument(
    "--template", type=int, required=True, help="default 0, template (error 1 in 4): 1 "
)
parser.add_argument(
    "--Lm", type=int, default=300, help="min no. snps to be considerred as IBD"
)
parser.add_argument("--Lf", type=float, default=3.0, help="min IBD length in cM")
parser.add_argument("--mem_gb", type=int, required=True)  # not used
parser.add_argument("--nthreads", type=int, default=None)
parser.add_argument("--minmac", type=int, default=1, help="min MAC")
parser.add_argument("--use_phase_correction", type=int, default=1, choices=[0, 1])
parser.add_argument("--genome_set_id", type=int, required=True)

args = parser.parse_args()

vcf = args.vcf
bp_per_cm = int(0.01 / args.r)
seqlen_in_cm = args.seqlen / bp_per_cm
chrno = args.chrno
mem_gb = args.mem_gb
nthreads = args.nthreads
template_opt = args.template
Lm = args.Lm
Lf = args.Lf
minmac = args.minmac
use_phase_correction = True if args.use_phase_correction == 1 else False

# filter by minmac and
# decompression vcf (required by phasedibd)
nsam = len(allel.read_vcf_headers(vcf).samples)
minmaf = minmac / nsam / 2
subprocess.run(
    f"bcftools view -q {minmaf}:minor {vcf} > tmp.vcf", shell=True, check=True
)
vcf = f"tmp.vcf"

# make a map for every snp position of the vcf file
data = allel.read_vcf(vcf, fields=["variants/POS", "variants/CHROM", "samples"])
samples = pd.Series(data["samples"])
pos = data["variants/POS"]
chrom = data["variants/CHROM"]
centimorgan = pos / bp_per_cm
gmap = pd.DataFrame({"Chrom": chrom, "Id": ".", "cM": centimorgan, "Bp": pos})
gmap.to_csv(f"{chrno}.map", sep=" ", index=None, header=None)

# haplotypes
haplotypes = phasedibd.VcfHaplotypeAlignment(vcf, f"{chrno}.map")

# define template
if template_opt == 0:
    tpbwt = phasedibd.TPBWTAnalysis()
elif template_opt == 1:
    template = [
        [0, 1, 1, 1],
        [1, 0, 1, 1],
        [1, 1, 0, 1],
        [1, 1, 1, 0],
    ]
    tpbwt = phasedibd.TPBWTAnalysis(template)
else:
    raise Exception("unknown template")

# call ibd
ibd_results = tpbwt.compute_ibd(
    haplotypes, L_f=Lf, L_m=Lm, use_phase_correction=use_phase_correction
)

# check ids
# tpbwt is using index not sample names
# check the vcf file names are in order so corrisponding to tpbwt ids
assert (
    samples.str.replace("tsk_", "", regex=False).astype(int) == np.arange(samples.size)
).all()

# remove duplicated ibd
ibd_results = ibd_results[lambda x: (x.id1_haplotype == 0) & (x.id2_haplotype == 0)]

# reformat
# Id1", "Id2", "Start", "End" , "Ancestor", "Tmrca",  "HasMutation"
ibd_df = pd.DataFrame(
    {
        "Id1": ibd_results["id1"],
        "Id2": ibd_results["id2"],
        "Start": ibd_results["start_bp"],
        "End": ibd_results["end_bp"],
        "Ancestor": 99999,
        "Tmrca": 100,
        "HasMutation": 0,
    }
)

# write to file
ofn = f"{args.genome_set_id}_{chrno}_tpbwtibd.ibd"
ibd_df.to_csv(ofn, sep="\t", header=True, index=None)

print(
    f"""
    output: 
        {ofn}
     """
)
# remove the temporiy vcf file
Path("tmp.vcf").unlink()
