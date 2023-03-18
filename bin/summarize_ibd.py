#! /usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import tskit
import pyranges as pr

parser = argparse.ArgumentParser()
parser.add_argument("--ibd_tsv", type=str, required=True)
if sys.argv[0] == "":
    arg_list = [
        "--ibd_tsv",
        "/autofs/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd/runs/run_test/results/N0_10000_Neutral/s7_summarize_ibd/info.tsv",
    ]
    args = parser.parse_args(arg_list)
else:
    args = parser.parse_args()

print(args)

ibd_tsv = args.ibd_tsv

tsinfer_node_time_thres = [300, 1000, 3000, 10000]


def get_ibd_coverage(ibd_fn, chrno: int, srclabel: str):
    # test
    # chrno = 1
    # ibd_fn = "/autofs/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd/runs/run_test/work/5f/2a3673d0ebceb2f241c9415dfb8d8a/1.ibd"
    df = pd.read_csv(ibd_fn, sep="\t")
    max_bp = df.End.max()
    num_points = 200
    sampling_points = np.linspace(100, max_bp - 100, num=num_points, dtype=int)
    sampling_points = pr.PyRanges(
        chromosomes=[chrno] * num_points,
        starts=sampling_points,
        ends=sampling_points + 1,
    )
    d = {}
    if srclabel != "tsinfer_tree":
        d[srclabel] = pr.PyRanges(
            chromosomes=np.repeat(chrno, df.shape[0]), starts=df.Start, ends=df.End
        )
    else:
        for thres in tsinfer_node_time_thres:
            df_temp = df[df.Tmrca <= thres]
            d[f"{srclabel}_upto_n{thres}"] = pr.PyRanges(
                chromosomes=np.repeat(chrno, df_temp.shape[0]),
                starts=df_temp.Start,
                ends=df_temp.End,
            )
    coverage = pr.count_overlaps(d, sampling_points)
    return coverage


info = pd.read_csv(ibd_tsv, sep="\t", names=["SimLabel", "SrcLabel", "ChrNo", "Path"])
info.sort_values(["SimLabel", "ChrNo", "SrcLabel"], inplace=True)
info

for simlabel in info.SimLabel.unique():
    png_fn = f"{simlabel}.png"
    sub_info = info[info.SimLabel == simlabel]
    if sub_info.shape[0] == 0:
        continue
    uniq_chnos = sub_info.ChrNo.unique()
    uniq_src_labels = sub_info.SrcLabel.unique()
    sz_chrno = uniq_chnos.size
    sz_src = uniq_src_labels.size
    num_cols = sz_src - 1 + len(tsinfer_node_time_thres)
    fig, axes = plt.subplots(
        nrows=sz_chrno,
        ncols=num_cols,
        figsize=(3 * num_cols, 3 * sz_chrno),
        sharey="col",
        constrained_layout=True,
    )
    fig.suptitle(simlabel)
    for i, chrno in enumerate(uniq_chnos):
        axes_row = axes if sz_chrno == 1 else axes[i]
        x_list = []
        y_list = []
        label_list = []
        for j, srclabel in enumerate(uniq_src_labels):
            sel_info = sub_info[
                (sub_info.ChrNo == chrno) & (sub_info.SrcLabel == srclabel)
            ]
            if sel_info.shape[0] == 0:
                continue
            ibd_fn = sel_info.Path.values[0]
            print(ibd_fn)
            coverage = get_ibd_coverage(ibd_fn, chrno, srclabel)
            for col in coverage.columns[3:]:
                label_list.append(col)
                x_list.append(coverage.Start)
                y_list.append(coverage.df[col])
        for i, (ax, label, x, y) in enumerate(
            zip(axes_row, label_list, x_list, y_list)
        ):
            ax.plot(x, y, "--.", label=label)
            ax.set_ylabel("IBD Coverage")
            ax.set_xlabel("Position in BP")
            ax.legend()
            ax.set_title(f"{chrno} | {label}")
    plt.savefig(png_fn)
    plt.close()
