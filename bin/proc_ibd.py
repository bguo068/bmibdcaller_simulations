#! /usr/bin/env python3
import numpy as np
import pandas as pd
import pyranges as pr
import pathlib
import argparse
import matplotlib.pyplot as plt
import ibdpeak

# ---------------------------------------------------------------------

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--chrno", type=int, default=1, help="chr number")
parser.add_argument("--cut_start", type=int, default=0, help="start pos in bp to cut")
parser.add_argument("--cut_end", type=int, default=1, help="end pos in bp to cut")
parser.add_argument("--bp_per_cm", type=int, default=1000_000, help="bp per cm")
parser.add_argument("--peaks_to_rm", type=str, default=None, help="peaks_to_rm")

parser.add_argument("--ibd_in", type=str, required=True, help="input ibd file")
choices = ["nocut", "autocut", "hardcut", "mutcut", 'autocut_hapibd']
parser.add_argument(
    "--ibd_cut_mode", type=str, choices=choices, default="nocut", help="cut mode"
)
args = parser.parse_args()
# args = parser.parse_args("--ibd_in 1.ibd".split())
print(args)

chrno = args.chrno
cut_start = args.cut_start
cut_end = args.cut_end
bp_per_cm = args.bp_per_cm
ibd_fn = args.ibd_in

# ---------------------------------------------------------------------
# this works for ibd file with single or multiple chromosomes
ibd_path = pathlib.Path(ibd_fn)
ibd_out = ibd_path.with_suffix(".ibd_proc")
map_out = ibd_path.with_suffix(".map_proc")
targets_out = ibd_path.with_suffix(".targets_cut")


# steps:
# 1. read raw ibd
# 2. make diploid
# 3. filter raw ibd like tmrca, whether contains mutations
# 4. remove homozygosity by decent
# 5. flatten ibd
# 6. cut ibd either automatically or by hard setting start and end points

# 1. loaded uncut IBD file
df_ibd = pd.read_csv(ibd_fn, sep="\t")

# 2.1 make diploid
df_ibd["Chromosome"] = chrno
df_ibd["Sample1"] = df_ibd["Id1"] // 2
df_ibd["Sample2"] = df_ibd["Id2"] // 2
df_ibd["Hap1"] = df_ibd["Id1"] % 2
df_ibd["Hap2"] = df_ibd["Id2"] % 2

# 2.2 ordering samples per record
sel1 = df_ibd.Sample1 > df_ibd.Sample2
sel2 = (df_ibd.Sample1 == df_ibd.Sample2) & (df_ibd.Hap1 < df_ibd.Hap2)
sel = sel1 | sel2
old_cols = ["Sample1", "Hap1", "Sample2", "Hap2"]
new_cols = ["Sample2", "Hap2", "Sample1", "Hap1"]
df_ibd.loc[sel, old_cols] = df_ibd.loc[sel, new_cols].to_numpy()

# 3-4. sel by trmca and HasMutation
sel = df_ibd.Tmrca > 1.5
if args.ibd_cut_mode == "mutcut":
    sel = sel & (df_ibd.HasMutation == 0)
sel = sel & (df_ibd.Sample1 != df_ibd.Sample2)

sel_cols = ["Sample1", "Sample2", "Chromosome", "Start", "End"]
df_ibd = df_ibd.loc[sel, sel_cols].copy()

# 5. merge ibd segment per sample pair
fake_chrom = df_ibd[["Sample1", "Sample2", "Chromosome"]].astype(str)
fake_chrom = fake_chrom["Sample1"].str.cat(
    fake_chrom[["Sample2", "Chromosome"]], sep=":"
)
gr = pr.PyRanges(chromosomes=fake_chrom, starts=df_ibd.Start, ends=df_ibd.End)
gr = gr.merge()

# reformat to normal IBD format
df_ibd_merged = gr.df[["Start", "End"]]
df_ibd_merged[["Sample1", "Sample2", "Chromosome"]] = gr.df.Chromosome.str.split(
    ":", expand=True
).astype(int)
gr = pr.PyRanges(df_ibd_merged)


# 6. cut-ibd
if args.ibd_cut_mode == "autocut":
    finder = ibdpeak.IbdPeaks(bp_per_cm, sampling_step=0.01 * bp_per_cm)
    target_regions, cores, cov = finder.find_peaks(gr)
elif args.ibd_cut_mode == "autocut_hapibd":
    tmp = pd.read_csv(args.peaks_to_rm, sep="\t", names=["Chromosome", "Start", "End"])
    target_regions = pr.PyRanges(tmp)
elif args.ibd_cut_mode == "hardcut":
    # make a range object containing the peak regions
    chrs = gr.df.Chromosome.unique().tolist()
    target_regions = pr.PyRanges(
        chromosomes=chrs, starts=[cut_start] * len(chrs), ends=[cut_end] * len(chrs)
    )
else:
    target_regions = pr.PyRanges(chromosomes=[], starts=[], ends=[])

# write target regions to a file
target_regions.df.to_csv(targets_out, sep="\t", index=None)
print(targets_out)

# make a complement ranges of the exclusion regions
all_chr_regions = pr.PyRanges(chromosomes=[chrno], starts=[0], ends=[gr.End.max()])
include_regions = all_chr_regions.subtract(target_regions)

# make each included region as a separate chromosome
gr_cut_list = []
df_map_list = []
for i, chrom, start, end in include_regions.df[
    ["Chromosome", "Start", "End"]
].itertuples():
    new_chrom = f"{chrom}_{i}"
    # ibd: make new chromosome for each included ranges
    single_range = pr.PyRanges(chromosomes=[chrom], starts=[start], ends=[end])
    gr_cut_single_range = gr.intersect(single_range)
    tmp_df = gr_cut_single_range.df.copy()
    tmp_df["Chromosome"] = new_chrom
    gr_cut_list.append(pr.PyRanges(tmp_df))

    # map
    df_map_single = pd.DataFrame(
        {
            "chromosome": new_chrom,
            "RID": ".",
            "Cm": [start / bp_per_cm, end / bp_per_cm],
            "Bp": [start, end],
        }
    )
    df_map_list.append(df_map_single)


# subtract
gr_cut = pr.concat(gr_cut_list)

# update cM fields (grange only update chromosome/start/end columns)

# add hap columns
df = gr_cut.df.copy()
df["Cm"] = (df.End - df.Start) / bp_per_cm
df["Hap1"] = 0
df["Hap2"] = 0
df = df[["Sample1", "Hap1", "Sample2", "Hap2"] + ["Chromosome", "Start", "End", "Cm"]]
df.to_csv(ibd_out, sep="\t", header=None, index=None)

print(ibd_out)

# output map
df_map = pd.concat(df_map_list, axis=0)
df_map.to_csv(map_out, sep=" ", index=None, header=None)

"""
# map already generated by tskibd
# df_map = trueibd.make_plink_map(chrno, bp_per_cm, seqlen=seqlen)
# df_map.to_csv(f"{chrno}.map", sep=" ", header=None, index=None)
"""
"""
output:

"{chrno}.ibd_proc"
"{chrno}.map_proc"
"{chrno}.targets_cut")
"""
