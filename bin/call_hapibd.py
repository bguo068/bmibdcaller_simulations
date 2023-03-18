#! /usr/bin/env python3
import argparse
from pathlib import Path
import sys
import tskit
import tsinfer
import tsdate
import msprime
import pandas as pd
import numpy as np
import allel
import gzip
from typing import Tuple, List
import io
import msprime
import sys
import subprocess

# --------------------- get hapibd jar path------------------------
if "__file__" in locals():
    hapibd_jar = str((Path(__file__).parents[1] / "lib/hap-ibd.jar").absolute())
else:
    hapibd_jar = "/autofs/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd/trueibd/nf_scripts/lib/hap-ibd.jar"

# --------------------- parse arguments -----------------------------
parser = argparse.ArgumentParser(
    """
call_hapibd.py does to three things
1. automatically generate the genetic map
2. call hapibd on a single-chromosome vcf file
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
parser.add_argument("--minseed", type=float, default=2.0)
parser.add_argument("--minoutput", type=float, default=2.0)
parser.add_argument("--maxgap", type=int, default=1000)
parser.add_argument("--minextend", type=float, default=1.0)
parser.add_argument("--minmarkers", type=int, default=100)
parser.add_argument("--mem_gb", type=int, required=True)
parser.add_argument("--nthreads", type=int, default=None)

if sys.argv[0] == "":  # interactive mode, used for testing
    args_lst = "--vcf 1.vcf.gz --bp_per_cm 15000 --seqlen_in_cm 100 --chrno 1 --mem_gb 10 --nthreads 10".split()
    args = parser.parse_args(args_lst)
else:
    args = parser.parse_args()

vcf = args.vcf
bp_per_cm = args.bp_per_cm
seqlen_in_cm = args.seqlen_in_cm
chrno = args.chrno
mem_gb = args.mem_gb
nthreads = args.nthreads
minseed = args.minseed
minoutput = args.minoutput
maxgap = args.maxgap
minextend = args.minextend
minmarkers = args.minmarkers

# ------------------- make genetic map ----------------------------
with open(f"{chrno}.map", "w") as f:
    fields = [f"{chrno}", ".", "0", "1"]
    f.write("\t".join(fields) + "\n")
    fields = [f"{chrno}", ".", f"{seqlen_in_cm}", f"{bp_per_cm * seqlen_in_cm}"]
    f.write("\t".join(fields) + "\n")

# ------------------- call hapibd ----------------------------------
map_fn = f"{chrno}.map"
cmd = (
    f"java -Xmx{mem_gb}g -jar {hapibd_jar} gt={vcf} map={map_fn} nthreads={nthreads} out={chrno}"
    f" min-seed={minseed} min-output={minoutput} max-gap={maxgap} min-extend={minextend} min-markers={minmarkers}"
)
print(cmd)
res = subprocess.run(cmd.split(), text=True, capture_output=True)
# check err
if res.returncode != 0 or "ERROR" in res.stderr.capitalize():
    raise Exception(res.stderr)

# ----------------- reformat hapibd results -----------------------
df_ibd = pd.read_csv(f"{chrno}.ibd.gz", sep="\t", header=None)
df_ibd.columns = ["Id1", "Hap1", "Id2", "Hap2", "Chrom", "Start", "End", "Cm"]
# remove replicative records and use needed columns
df_ibd = df_ibd.loc[
    lambda df: (df.Hap1 == 1) & (df.Hap2 == 1), ["Id1", "Id2", "Start", "End"]
]
# update tsk sample name to nunber only names
df_ibd["Id1"] = df_ibd["Id1"].str.replace("tsk_", "", regex=False)
df_ibd["Id2"] = df_ibd["Id2"].str.replace("tsk_", "", regex=False)
# add fake columns
df_ibd["Ancestor"] = 99999
df_ibd["Tmrca"] = 100
df_ibd["HasMutation"] = 0

# ---------------- output IBD -------------------------------
df_ibd.to_csv(f"{chrno}.ibd", header=True, index=None, sep="\t")
