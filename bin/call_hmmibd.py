#! /usr/bin/env python3
import allel
from subprocess import run
from pathlib import Path
import pandas as pd
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=str, required=True)
    help = (
        "Note: bp_per_cm is not respected by hmmIBD. Recombination rate is"
        "hard-coded in hmmIBD.c source code"
    )
    parser.add_argument("--r", type=float, required=True, help=help)
    parser.add_argument("--seqlen", type=int, required=True)
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--chrno", type=int, required=True)
    parser.add_argument("--minmac", type=int, default=1, help="min MAC")
    parser.add_argument("--mincm", type=float, default=2.0)
    parser.add_argument("--genome_set_id", type=int, required=True)

    return parser.parse_args()


args = get_args()
vcf = args.vcf
bp_per_cm = 0.01 / args.r
seqlen_in_cm = args.seqlen / bp_per_cm
chrno = args.chrno
n = args.n
m = args.m
minmac = args.minmac
mincm = args.mincm

# ------------------- make genetic map ----------------------------
with open(f"{chrno}.map", "w") as f:
    fields = [f"{chrno}", ".", "0", "1"]
    f.write("\t".join(fields) + "\n")
    fields = [f"{chrno}", ".", f"{seqlen_in_cm}", f"{bp_per_cm * seqlen_in_cm}"]
    f.write("\t".join(fields) + "\n")

# ------------------- prepare hmmibd input ----------------------------
calldata = allel.read_vcf(vcf, fields=["samples", "CHROM", "POS", "GT"])
samples = calldata["samples"]
gt = calldata["calldata/GT"][:, :, 0]
minmac = minmac // 2  # given GT is using ploidy 0 only in the above line

df_pos = pd.DataFrame(
    {
        "chrom": calldata["variants/CHROM"],
        "pos": calldata["variants/POS"],
    }
)
df_gt = pd.DataFrame(gt, columns=samples)

# filter site by minmac
nsample = df_gt.shape[1]
ac = df_gt.sum(axis=1).to_numpy()
# .  get mac from ac
sel = ac * 2 > nsample
ac[sel] = nsample - ac[sel]
sel = ac >= minmac
df_pos = df_pos.loc[sel, :]
df_gt = df_gt.loc[sel, :]


df_hmm_inputs = pd.concat([df_pos, df_gt], axis=1)
df_hmm_inputs.to_csv("hmm_inputs.txt", sep="\t", index=None)

# ------------------- run hmmibd  ----------------------------
# cmd = f"hmmIBDr -i hmm_inputs.txt -o hmmibd_out -n {n} -m {m} -r {args.r}"
cmd = f"hmmibd2 -i hmm_inputs.txt -o hmmibd_out -n {n} -m {m} -r {args.r}"
run(cmd, shell=True)


df_ibd = pd.read_csv("hmmibd_out.hmm.txt", sep="\t")[lambda x: (x["different"] == 0)]
df_ibd.columns = ["Id1", "Id2", "Chr", "Start", "End", "Diff", "Nsnp"]
# update tsk sample name to nunber only names
df_ibd["Id1"] = df_ibd["Id1"].str.replace("tsk_", "", regex=False)
df_ibd["Id2"] = df_ibd["Id2"].str.replace("tsk_", "", regex=False)
# add fake columns
df_ibd["Ancestor"] = 99999
df_ibd["Tmrca"] = 100
df_ibd["HasMutation"] = 0

sel_cols = ["Id1", "Id2", "Start", "End", "Ancestor", "Tmrca", "HasMutation"]
sel_rows = ((df_ibd.End - df_ibd.Start) / bp_per_cm) >= mincm

ofn = f"{args.genome_set_id}_{chrno}_hmmibd.ibd"
df_ibd.loc[sel_rows, sel_cols].to_csv(ofn, sep="\t", index=None)

# -------------------- write log (just for consistency with other ibd callers) -
Path(f"{chrno}.log").write_text(cmd)

print(
    f"""
    output: 
        {ofn}
     """
)
