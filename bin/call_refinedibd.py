#! /usr/bin/env python3
import argparse
import subprocess
from pathlib import Path
import allel

import pandas as pd

# --------------------- get hapibd jar path------------------------
refinedibd_jar = str(
    (Path(__file__).parents[1] / "lib/refined-ibd.17Jan20.102.jar").resolve()
)

# --------------------- parse arguments -----------------------------
parser = argparse.ArgumentParser()

parser.add_argument("--vcf", type=str, required=True, help="input vcf file")
parser.add_argument("--r", type=float, required=True)
parser.add_argument("--seqlen", type=int, required=True)
parser.add_argument("--chrno", type=int, required=True)
parser.add_argument("--window", type=float, default=40.0)
parser.add_argument("--lod", type=float, default=0.3)
parser.add_argument("--length", type=float, default=2.0)
parser.add_argument("--scale", type=float, default=0)
parser.add_argument("--mem_gb", type=int, required=True)
parser.add_argument("--minmac", type=int, default=1, help="min MAC")
parser.add_argument("--nthreads", type=int, default=None)
parser.add_argument("--genome_set_id", type=int, required=True)

args = parser.parse_args()

vcf = args.vcf
bp_per_cm = int(0.01 / args.r)
seqlen_in_cm = args.seqlen / bp_per_cm
chrno = args.chrno
window = args.window
lod = args.lod
length = args.length
scale = args.scale
minmac = args.minmac
mem_gb = args.mem_gb
nthreads = args.nthreads

# filter by minmac and
# decompression vcf (required by phasedibd)
nsam = len(allel.read_vcf_headers(vcf).samples)
minmaf = minmac / nsam / 2
subprocess.run(
    f"bcftools view -q {minmaf}:minor {vcf} -Oz -o tmp.vcf.gz", shell=True, check=True
)
vcf = f"tmp.vcf.gz"


# ------------------- make genetic map ----------------------------
with open(f"{chrno}.map", "w") as f:
    fields = [f"{chrno}", ".", "0", "1"]
    f.write("\t".join(fields) + "\n")
    fields = [f"{chrno}", ".", f"{seqlen_in_cm}", f"{int(bp_per_cm * seqlen_in_cm)}"]
    f.write("\t".join(fields) + "\n")

# ------------------- call hapibd ----------------------------------
map_fn = f"{chrno}.map"
cmd = (
    f"java -Xss5M -Xmx{mem_gb}g -jar {refinedibd_jar} gt={vcf} map={map_fn}"
    f" nthreads={nthreads} out={chrno} chrom={chrno}"
    f" lod={lod} window={window} scale={scale} length={length}"
)
print(cmd)
res = subprocess.run(cmd.split(), text=True, capture_output=True)
# check err
if (
    res.returncode != 0
    or "ERROR" in res.stderr.capitalize()
    or "ERROR" in res.stdout.capitalize()
):
    raise Exception(res.stderr + "\n" + res.stdout)

# ----------------- reformat refinedibd results (there are 9 columns) ---------------
df_ibd = pd.read_csv(f"{chrno}.ibd.gz", sep="\t", header=None)
df_ibd.columns = ["Id1", "Hap1", "Id2", "Hap2", "Chrom", "Start", "End", "LOD", "Cm"]
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
ofn = f"{args.genome_set_id}_{chrno}_refinedibd.ibd"
df_ibd.to_csv(ofn, header=True, index=None, sep="\t")

print(
    f"""
    output: 
        {ofn}
     """
)

Path("tmp.vcf.gz").unlink()
