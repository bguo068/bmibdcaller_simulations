#! /usr/bin/env python3

"""
I had a lot of problems with getting quality IBD from sequencing data for
calling IBDNe.  I used to compare different IBD calling tools or parameters by
comparing the final Ne estimates.  It was straightforward but also very slow. 

The aim of this script is try to use shortcuts --  do quality checks on IBD
intead of directly looking at Ne.

1. one is to calculate false postive rate and false negative rate (hap-ibd
paper, Zhou et al 2020)
2. second is to mimicking strategy used in IBDNe-- binning the observed IBD
segments by their length and number of ends reaching the end of the chromosome.
(IBDNe paper, Browning et al 2015, p415)
3. three is to add error in total IBD number and length (see TPBWT paper,
Freyman et. al 2020)

"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path
from subprocess import run
import sys
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

parser = argparse.ArgumentParser()
parser.add_argument("--raw_true_ibd", type=str, nargs="+", required=True)
parser.add_argument("--raw_detected_ibd", type=str, nargs="+", required=True)
parser.add_argument("--bp_per_cm", type=int, required=True)
parser.add_argument("--genome_size_cm", type=float, required=True)
parser.add_argument("--out_prefix", type=str, required=True)
parser.add_argument("--min_tmrca", type=int, default=-1)
parser.add_argument("--max_tmrca", type=int, default=1_000_000_000)

if sys.argv[0] == "":  # for debug
    args_list = [
        "--raw_true_ibd",
        "/autofs/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd/"
        "runs/run_real/results/N0_10000_selCoeff_0.4_selStartG_50/s3_calc_raw_ibd/true_tree/1/1.ibd",
        "/autofs/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd/"
        "runs/run_real/results/N0_10000_selCoeff_0.4_selStartG_50/s3_calc_raw_ibd/true_tree/2/2.ibd",
        "--raw_detected_ibd",
        "/autofs/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd/"
        "runs/run_real/results/N0_10000_selCoeff_0.4_selStartG_50/s3_seqbased_raw_ibd/hapibd/1/1.ibd",
        "/autofs/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd/"
        "runs/run_real/results/N0_10000_selCoeff_0.4_selStartG_50/s3_seqbased_raw_ibd/hapibd/2/2.ibd",
        # "runs/run_real/results/N0_10000_selCoeff_0.4_selStartG_50/s3_seqbased_raw_ibd/refinedibd/1/1.ibd",
        # "runs/run_real/results/N0_10000_selCoeff_0.4_selStartG_50/s3_calc_raw_ibd/tsinfer_tree/1/1.ibd",
        "--bp_per_cm",
        "15000",
        "--genome_size_cm",
        "200",
        "--out_prefix",
        "./",
    ]
    args = parser.parse_args(args_list)
else:
    args = parser.parse_args()

raw_true_ibd = args.raw_true_ibd
raw_detected_ibd = args.raw_detected_ibd
bp_per_cm = args.bp_per_cm
out_prefix = args.out_prefix
genome_size_cm = args.genome_size_cm
min_tmrca = args.min_tmrca
max_tmrca = args.max_tmrca

# ------------------------ Get False Pos/Neg Rate -----------------------------------------------

# pyranges intersect can do the work but is much slower than bedtools
# reformat the files and then run bedtools intersect

n_files = len(raw_true_ibd)
assert n_files == len(raw_detected_ibd)

for i in range(n_files):
    cmd = f"""
    sed 1d {raw_true_ibd[i]} | \\
        awk -v OFS='\t' '$6 >= {min_tmrca} && $6 <= {max_tmrca} {{if($1>$2) {{id1=$1; id2=$2}} else {{id1=$2; id2=$1}}; print id1 ":" id2, $3, $4}}' | \\
        sort -k1,1 -k2,2n | \\
        bedtools merge -i - | \\
        awk -v OFS='\t' '{{print $1,$2,$3,($3-$2)/{bp_per_cm} }}' > true_{i}.bed
    sed 1d {raw_detected_ibd[i]} | \\
        awk -v OFS='\t' '{{if($1>$2) {{id1=$1; id2=$2}} else {{id1=$2; id2=$1}}; print id1 ":" id2, $3, $4}}' | \\
        sort -k1,1 -k2,2n | \\
        bedtools merge -i - | \\
        awk -v OFS='\t' '{{print $1,$2,$3,($3-$2)/{bp_per_cm} }}' > infer_{i}.bed
    bedtools intersect -a true_{i}.bed -b infer_{i}.bed > true_intersect_infer_{i}.bed
    bedtools intersect -a infer_{i}.bed -b true_{i}.bed > infer_intersect_true_{i}.bed
    """
    run(cmd, shell=True)


# red all bed files
def read_bed_file(prefix, n):
    return pd.concat(
        pd.read_csv(
            f"{prefix}_{i}.bed", sep="\t", names=["Chromosome", "Start", "End", "Cm"]
        )  # .assign(Chromosome=lambda x: x.Chromosome + ":" + str(i))
        for i in range(n)
    )


df_true = read_bed_file("true", n_files)
df_infer = read_bed_file("infer", n_files)
df_true_intersect_infer = read_bed_file("true_intersect_infer", n_files)
df_infer_intersect_true = read_bed_file("infer_intersect_true", n_files)

# assign bin by looking at original segment length
bins = [2.5, 3, 4, 6, 10, 18, np.inf]

for df in [df_true, df_infer, df_true_intersect_infer, df_infer_intersect_true]:
    df["Bin"] = pd.cut(df.Cm, bins, right=False)

# update intersection segment Cm
df_true_intersect_infer["Cm"] = (
    df_true_intersect_infer.End - df_true_intersect_infer.Start
) / bp_per_cm
df_infer_intersect_true["Cm"] = (
    df_infer_intersect_true.End - df_infer_intersect_true.Start
) / bp_per_cm

# false positive
#
df_false_pos = pd.concat(
    [
        df_infer_intersect_true.groupby("Bin")["Cm"].sum().rename("InferIntersectTrue"),
        df_infer.groupby("Bin")["Cm"].sum().rename("Infer"),
    ],
    axis=1,
).assign(
    FalsePosRate=lambda x: (x.Infer - x.InferIntersectTrue) / x.Infer.replace(0, np.nan)
)

# False negative
df_false_neg = pd.concat(
    [
        df_true_intersect_infer.groupby("Bin")["Cm"].sum().rename("TrueIntersectInfer"),
        df_true.groupby("Bin")["Cm"].sum().rename("True"),
    ],
    axis=1,
).assign(
    FalseNegRate=lambda x: (x["True"] - x.TrueIntersectInfer)
    / x["True"].replace(0, np.nan)
)

df_false_pos.to_csv(out_prefix + "false_pos.tsv", sep="\t", index=None)
df_false_neg.to_csv(out_prefix + "false_neg.tsv", sep="\t", index=None)

# ------------------------ Get totalIBD error -------------------------------------
# see TPBWT paper
# Note: here is chromosome names are just catenation of two sample id2
# so groupby generate genome wide aggregate when there are multiple real chromosomes

genome_wide_error = (
    pd.concat(
        [
            df_true.groupby("Chromosome")["Cm"]
            .agg(["count", sum])
            .rename(columns={"count": "NumTrueIBD", "sum": "TrueTotalIBD"}),
            df_infer.groupby("Chromosome")["Cm"]
            .agg(["count", sum])
            .rename(columns={"count": "NumInferIBD", "sum": "InferTotalIBD"}),
        ],
        axis=1,
    )
    .fillna(0)
    .assign(
        NumIBDError=lambda x: x.NumInferIBD - x.NumTrueIBD,
        TotalIBDError=lambda x: (x.InferTotalIBD - x.TrueTotalIBD) / genome_size_cm,
    )
)

genome_wide_error.to_csv(out_prefix + "chr_wide_error.tsv", sep="\t")
genome_wide_error.NumIBDError.describe()
genome_wide_error.TotalIBDError.describe()

# ------------------------ Get IBDNe estimation -----------------------------------------
win_size_cm = 0.05
df_true = read_bed_file("true", n_files)
df_infer = read_bed_file("infer", n_files)

max_bp = max(df_true.End.max(), df_infer.End.max()) - 0.2 * bp_per_cm
min_bp = 0.2 * bp_per_cm

df_true["End"] = df_true["End"].where(df_true.End < max_bp, max_bp)
df_true["Start"] = df_true["Start"].where(df_true.Start > min_bp, min_bp)
df_true["NumEnds"] = (df_true.End == max_bp).astype(int) + (
    df_true.Start == min_bp
).astype(int)
bins = pd.cut(
    df_true.Cm, bins=np.arange(min_bp / bp_per_cm, max_bp / bp_per_cm, win_size_cm)
)
df_true["Bin"] = bins.cat.categories.mid.values[bins.cat.codes]

df_infer["End"] = df_infer["End"].where(df_infer.End < max_bp, max_bp)
df_infer["Start"] = df_infer["Start"].where(df_infer.Start > min_bp, min_bp)
df_infer["NumEnds"] = (df_infer.End == max_bp).astype(int) + (
    df_infer.Start == min_bp
).astype(int)
bins = pd.cut(
    df_infer.Cm, bins=np.arange(min_bp / bp_per_cm, max_bp / bp_per_cm, win_size_cm)
)
# need the values of bins.cat.codes.values instead of the pandas series itself
df_infer["Bin"] = bins.cat.categories.mid.values[bins.cat.codes.values]

df_bin_totalibd = (
    pd.concat(
        [
            df_true.groupby(["NumEnds", "Bin"])["Cm"]
            .count()
            .rename("TrueIntervalIBDCount"),
            df_infer.groupby(["NumEnds", "Bin"])["Cm"]
            .count()
            .rename("InferIntervalIBDCount"),
        ],
        axis=1,
    )
    .fillna(0)
    .reset_index()
    .assign(
        TrueIntervalTotalIBD=lambda x: x.Bin * x.TrueIntervalIBDCount,
        InferIntervalTotalIBD=lambda x: x.Bin * x.InferIntervalIBDCount,
    )
)

df_bin_totalibd.to_csv(out_prefix + "binned_totalibd.tsv", sep="\t", index=None)


# ----------------------- Plotting -----------------------------
fig, axes = plt.subplots(ncols=4, figsize=(16, 4), constrained_layout=True)
# fn and pn
ax = axes[0]
df_plot = pd.concat(
    [df_false_pos["FalsePosRate"], df_false_neg["FalseNegRate"]], axis=1
)
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.plot(df_plot.FalsePosRate, df_plot.FalseNegRate)
for Bin, fp, fn in df_plot.itertuples():
    Bin = str(Bin)
    ax.plot([fp], [fn], ".", markersize=10, label=Bin)
ax.legend()
ax.set_ylabel("False-negative rate")
ax.set_xlabel("False-positive rate")
ax.set_title("IBD detection accuracy")
#  Interval Total IBD Count
ax = axes[1]
df_plot = df_bin_totalibd[
    lambda x: (x.TrueIntervalIBDCount > 0) & (x.InferIntervalIBDCount)
]
if df_plot.shape[0] > 20:
    q1, q2 = np.quantile(df_plot.Bin, q=[0.025, 0.975])
    df_plot = df_plot[df_plot.Bin.between(q1, q2, "both")]
    x = df_plot.Bin.values
    y = np.log2(df_plot.InferIntervalIBDCount / df_plot.TrueIntervalIBDCount)
    if x.size > 10:
        reg = LinearRegression().fit(x.reshape(-1, 1), y)
        a, b = reg.coef_, reg.intercept_
        y_pred = a * x + b
        ax.plot(x, y_pred, "r--", linewidth=2, label="observed (linear regression)")
    ax.plot(x, y, ".", label="log2 Inferred to True ratio")
    ax.plot(x, y != y, "k--", linewidth=2, label="expected")
    ax.legend()
ax.set_xlabel("Bin Center (cM)")
ax.set_ylabel("log2(Inferred Count/True Count")
ax.set_title(f"IBD Counts per length bin (width={win_size_cm}cM)")
# chromosome-wise error
ax = axes[2]
xx = genome_wide_error.TotalIBDError
if xx.size > 20:
    q1, q2 = np.quantile(xx, q=[0.025, 0.975])
    xx = xx[lambda x: x.between(q1, q2, "both")]
ax.hist(xx, bins=50)
ax.set_title("Total IBD Error / Chrom size")
ax.set_xlabel("Total IBD Error / chrom size")
ax.set_ylabel("Count")
# IBD number error
ax = axes[3]
q1, q2 = np.quantile(genome_wide_error.NumIBDError, q=[0.025, 0.975])
x, y = np.unique(
    genome_wide_error.NumIBDError[lambda x: x.between(q1, q2, "both")],
    return_counts=True,
)
ax.bar(x, y)
ax.set_title("Error in no. IBD per sample pair")
ax.set_xlabel("Error in no. IBD per sample pair")
ax.set_ylabel("Count")
plt.savefig(out_prefix + "plot.png")
