#! /usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("--ne_tsv", type=str, required=True)
parser.add_argument("--true_ne_tsv", type=str, required=True)
if sys.argv[0] == "":
    arg_list = [
        "--ne_tsv",
        "/data/bing/20220406_cmp_ibd_methods/run_test/work/0d/4810a70bdf88fcd034fceb41ae9759/info.tsv",
        "--true_ne_tsv",
        "/data/bing/20220406_cmp_ibd_methods/run_test/work/0d/4810a70bdf88fcd034fceb41ae9759/input.1",
    ]
    args = parser.parse_args(arg_list)
else:
    args = parser.parse_args()

ne_tsv = args.ne_tsv
true_ne_tsv = args.true_ne_tsv

info = pd.read_csv(
    ne_tsv, sep="\t", names=["SimLabel", "SrcLabel", "ProcLabel", "Path"]
)
info.sort_values(["SimLabel", "ProcLabel", "SrcLabel"], inplace=True)

truene = pd.read_csv(true_ne_tsv, sep="\t", names=["SimLabel", "TrueNePath"])

info = info.merge(truene, how="left", on="SimLabel")


for simlabel in info.SimLabel.unique():
    png_fn = f"{simlabel}.png"
    sub_info = info[info.SimLabel == simlabel]
    if sub_info.shape[0] == 0:
        continue
    uniq_proc_labels = sub_info.ProcLabel.unique()
    uniq_src_labels = sub_info.SrcLabel.unique()
    sz_proc = uniq_proc_labels.size
    sz_src = uniq_src_labels.size
    fig, axes = plt.subplots(
        nrows=sz_proc,
        ncols=sz_src,
        figsize=(3 * sz_src, 3 * sz_proc),
        constrained_layout=True,
    )
    fig.suptitle(simlabel)
    for i, proclabel in enumerate(uniq_proc_labels):
        for j, srclabel in enumerate(uniq_src_labels):
            if sz_proc == 1 and sz_src == 1:
                ax = axes
            if sz_proc == 1:
                ax = axes[j]
            elif sz_src == 1:
                ax = axes[i]
            else:
                ax = axes[i][j]
            sel_info = sub_info[
                (sub_info.ProcLabel == proclabel) & (sub_info.SrcLabel == srclabel)
            ]
            if sel_info.shape[0] == 0:
                continue
            ne_fn = sel_info.Path.values[0]
            truene_fn = sel_info.TrueNePath.values[0]
            # only read and plot when the file is not a place holder
            if Path(ne_fn).stat().st_size > 0:
                df = pd.read_csv(ne_fn, sep="\t")
                x = df.GEN
                y = df.NE
                y1 = df["LWR-95%CI"]
                y2 = df["UPR-95%CI"]
                ax.plot(x, y, "b", label="Estimated Ne")
                ax.fill_between(x, y1, y2, alpha=0.3)
            df_true = pd.read_csv(truene_fn, sep="\t")
            ax.plot(df_true.generation, df_true.true_ne, "r--", label="True Ne")
            ax.set_yscale("log")
            ax.set_ylim(10, 1e7)
            ax.set_xlim(0, 200)
            ax.set_xlabel("gen ago")
            ax.set_ylabel("Ne")
            ax.legend()
            ax.set_title(f"{proclabel}\n{srclabel}")
    plt.savefig(png_fn, dpi=600)
    plt.close()
