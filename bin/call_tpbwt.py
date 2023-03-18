#! /usr/bin/env python3

import phasedibd
import pathlib
import pandas as pd
import numpy as np
import allel
import argparse
from pathlib import Path
import sys
import pandas as pd
from typing import Tuple, List
import sys
import subprocess

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
parser.add_argument("--bp_per_cm", type=int, required=True)
parser.add_argument("--seqlen_in_cm", type=int, required=True)
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

if sys.argv[0] == "":  # interactive mode, used for testing
    args_lst = [
        "--vcf",
        "/local/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd/runs/run_0407/results/N0_10000_Neutral/s1_simulations/1/1.vcf.gz",
        "--template",
        "0",
        "--bp_per_cm",
        "15000",
        "--seqlen_in_cm",
        "100",
        "--chrno",
        "1",
        "--Lm",
        "300",
        "--Lf",
        "3.0",
        "--mem_gb",
        "10",
        "--nthreads",
        "10",
    ]
    args = parser.parse_args(args_lst)
else:
    args = parser.parse_args()

vcf = args.vcf
bp_per_cm = args.bp_per_cm
seqlen_in_cm = args.seqlen_in_cm
chrno = args.chrno
mem_gb = args.mem_gb
nthreads = args.nthreads
template_opt = args.template
Lm = args.Lm
Lf = args.Lf

# decompress vcf file if needed
if Path(vcf).suffix == ".gz":
    subprocess.run(f"gzip -dc {vcf} > {chrno}.vcf", shell=True)
    vcf = f"{chrno}.vcf"

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
ibd_results = tpbwt.compute_ibd(haplotypes, L_f=Lf, L_m=Lm)

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
out_ibd_fn = f"{chrno}.ibd"
ibd_df.to_csv(out_ibd_fn, sep="\t", header=True, index=None)

# write log file
Path(f"{chrno}.log").write_text(f"ibd no.: {ibd_df.shape[0]}\n")
