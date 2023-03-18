#! /usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import tskit

parser = argparse.ArgumentParser()
parser.add_argument("--ts_tsv", type=str, required=True)
if sys.argv[0] == "":
    arg_list = [
        "--ts_tsv",
        "/autofs/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd"
        + "/runs/run_test/results/N0_10000_Neutral/s8_summarize_treeseq/info.tsv",
    ]
    args = parser.parse_args(arg_list)
else:
    args = parser.parse_args()

print(args)

ts_tsv = args.ts_tsv


def get_tree_stats(ts: tskit.TreeSequence) -> dict:
    ## test
    # tree_fn = (
    #     "/autofs/projects-t3/toconnor_grp/bing.guo/temp_work/20220316_tsinfer_ibd_vs_hapibd"
    #     + "/runs/run_test/work/74/5cd0f1ade8a77de2e4a9429795cbe0/1.trees"
    # )
    # ts = tskit.load(tree_fn)

    print("in")
    tree_centers = []
    tree_spans = []
    tree_heights = []
    for tree in ts.trees():
        tree_centers.append((tree.interval.right + tree.interval.left) / 2)
        tree_spans.append(tree.span)
        tree_heights.append(tree.time(tree.root))
    tree_centers = np.array(tree_centers)
    tree_spans = np.array(tree_spans)
    tree_heights = np.array(tree_heights)

    step = int(ts.sequence_length / 20)
    df_tree_heights = pd.DataFrame(
        {
            "Center": tree_centers,
            "Span": tree_spans,
            "Height": tree_heights,
            "Window": np.round(tree_centers / step),
        }
    )
    df_tree_heights = (
        df_tree_heights.groupby("Window")
        .apply(lambda df: np.average(df.Height, weights=df.Span))
        .rename("AvgHeight")
        .reset_index()
    )
    df_tree_heights["WinCenter"] = (df_tree_heights["Window"] + 0.5) * step

    # avg_tree_span = tree_spans.mean()
    # weighted_tree_height = np.average(tree_heights, weights=tree_spans)

    stats = {
        "num_trees": ts.num_trees,
        # "avg_tree_span": avg_tree_span,
        "num_edges": ts.num_edges,
        "num_nodes": ts.num_nodes,
        # "avg_tree_height": weighted_tree_height,
        "node_times": ts.tables.nodes.time,
        "window_tree_heights": df_tree_heights,
    }
    print("out")

    return stats


info = pd.read_csv(ts_tsv, sep="\t", names=["SimLabel", "SrcLabel", "ChrNo", "Path"])
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
    num_stats = 5
    fig, axes = plt.subplots(
        nrows=sz_chrno,
        ncols=num_stats,
        figsize=(3 * num_stats, 3 * sz_chrno),
        constrained_layout=True,
    )
    fig.suptitle(simlabel)
    for i, chrno in enumerate(uniq_chnos):
        axes_row = axes if sz_chrno == 1 else axes[i]
        ss_list = []
        src_labels = []
        for j, srclabel in enumerate(uniq_src_labels):
            sel_info = sub_info[
                (sub_info.ChrNo == chrno) & (sub_info.SrcLabel == srclabel)
            ]
            if sel_info.shape[0] == 0:
                continue
            ts_fn = sel_info.Path.values[0]
            print(f"stats for {chrno}-{srclabel}...")
            ts = tskit.load(ts_fn)
            stats = get_tree_stats(ts)
            ss_list.append(stats)
            src_labels.append(srclabel)
        # --------------
        ylabel = "num_trees"
        ax = axes_row[0]
        y = [stats[ylabel] for stats in ss_list]
        ax.bar(src_labels, y)
        ax.set_ylabel(ylabel)
        ax.set_yscale("log")
        ax.set_ylim(1, 2 * max(y))
        ax.set_title(f"{chrno} | {ylabel}")
        # --------------
        ylabel = "num_edges"
        ax = axes_row[1]
        y = [stats[ylabel] for stats in ss_list]
        ax.bar(src_labels, y)
        ax.set_ylabel(ylabel)
        ax.set_yscale("log")
        ax.set_ylim(1, 2 * max(y))
        ax.set_title(f"{chrno} | {ylabel}")
        # --------------
        ylabel = "num_nodes"
        ax = axes_row[2]
        y = [stats[ylabel] for stats in ss_list]
        ax.bar(src_labels, y)
        ax.set_ylabel(ylabel)
        ax.set_yscale("log")
        ax.set_ylim(1, 2 * max(y))
        ax.set_title(f"{chrno} | {ylabel}")
        ax = axes_row[3]
        for srclabel, ss in zip(src_labels, ss_list):
            df = ss["window_tree_heights"]
            xlabel = "win_center"
            ylabel = "win_tree_height"
            x = df.WinCenter
            y = df.AvgHeight
            ax.plot(x, y, "--o", label=srclabel)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend()
            ax.set_title(f"{chrno} | win_tree_heights")
        # --------------
        ax = axes_row[4]
        for srclabel, ss in zip(src_labels, ss_list):
            name = "node_times"
            xlabel = "log10 ancestor node time"
            ylabel = "density"
            x = ss[name]
            x = np.log10(x[x > 0])
            ax.hist(x, histtype="step", density=True, bins=100, label=srclabel)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.legend()
            ax.set_title(f"{chrno} | node times")
    plt.savefig(png_fn)
    plt.close()
