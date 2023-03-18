#! /usr/bin/env python3

import pybedtools as pb
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--raw_ibd", type=str, required=True)
parser.add_argument("--chr_name", type=str, required=True)
parser.add_argument("--chr_start", type=int, default=0)
parser.add_argument("--step", type=int, default=150)
parser.add_argument("--chr_end", type=int, default=0)
parser.add_argument("--out_bed", type=str, default="out.bed")
args = parser.parse_args(
    # for testing
    # [
    #     "--raw_ibd",
    #     "/local/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd/runs/run0426_maxgap300_lf200_nemincm4/results/N0_10000_selCoeff_0.4_selStartG_50/s3_seqbased_raw_ibd/hapibd/1/1.ibd",
    #     # "/local/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd/runs/run0426_maxgap300_lf200_nemincm4/results/N0_10000_Neutral/s3_seqbased_raw_ibd/hapibd/1/1.ibd",
    #     "--chr_name",
    #     "1",
    #     "--chr_start",
    #     "1",
    #     "--chr_end",
    #     "0",
    #     "--step",
    #     "150",
    #     "--out_bed",
    #     "1_peaks.bed",
    # ]
)

ibd_fn = args.raw_ibd
chr_name = args.chr_name
chr_start = args.chr_start
chr_end = args.chr_end
step = args.step
out_bed = args.out_bed

# convert ibd data to data frame
df = pd.read_csv(ibd_fn, usecols=[2, 3], sep="\t")
df.insert(0, "chrom", chr_name)
ibd_bed = pb.BedTool.from_dataframe(df)
max_end = df["End"].max() + 1
if chr_end == 0:
    chr_end = max_end
del df

# make sampling points
x = np.arange(chr_start + step // 2, chr_end, step=step, dtype=int)
df = pd.DataFrame({"chrom": chr_name, "start": x, "end": x + 1})
sp_bed = pb.BedTool.from_dataframe(df)
del x, df


# coverage
cov_bed = sp_bed.intersect(ibd_bed, c=True)
cov_df = cov_bed.to_dataframe().rename(columns={"name": "Cov"})

# threshold
q5, q25, q50, q75, q95 = np.quantile(cov_df.Cov, q=[0.05, 0.25, 0.5, 0.75, 0.95])
iqr = q75 - q25
trim_mean = cov_df.Cov[lambda s: (s >= q5) & (s < q95)].mean()
trim_std = cov_df.Cov[lambda s: (s >= q5) & (s < q95)].std()

# core regions
core_df = cov_df.loc[
    lambda x: (x.Cov > q50 + 1.5 * iqr)  # & (x.Cov > trim_mean + 2 * trim_std)
]
core_bed = pb.BedTool.from_dataframe(core_df.iloc[:, :3]).merge(d=step)

# extension region
ext_df = cov_df.loc[lambda x: x.Cov > q50]
ext_bed = pb.BedTool.from_dataframe(ext_df.iloc[:, :3]).merge(d=step)

# extension region intersect with core regions
peaks = (
    ext_bed.intersect(core_bed, wa=True).merge(d=step).saveas(out_bed).to_dataframe()
)

# plot peaks and ibd coverage
_, ax = plt.subplots(constrained_layout=True)
ax.plot(cov_df["start"], cov_df.Cov, "r.")
ax.set_ylabel("IBD coverage")
ax.set_xlabel("Position(bp)")
for s, e in peaks.iloc[:, 1:].itertuples(index=False):
    sub = cov_df[(cov_df["start"] >= s) & (cov_df["end"] <= e)]
    ax.fill_between(sub["start"], y1=q50, y2=sub.Cov)
plt.savefig(Path(out_bed).with_suffix(".png"))
