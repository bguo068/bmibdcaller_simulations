#! /usr/bin/env python3
import tomli_w
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from subprocess import run
from shutil import copyfile, rmtree
import numpy as np
import argparse

p = argparse.ArgumentParser()
p.add_argument("id", type=int, choices=list(range(4)))
which = p.parse_args().id

# groups
demog = ["sp", "mp", "uknogc", "uknogcpf"][which]
r = 0.01 / ([15000] * 2 + [1000000] * 1 + [15000] * 1)[which]
len_bp = ([15000 * 100] * 2 + [1000000 * 60] * 1 + [15000 * 60] * 1)[which]
len_cm = ([100.0] * 2 + [60.0] * 1 + [60.0] * 1)[which]
genome_set_name = ["sp_neu", "mp_s00", "uk_gc0", "uk_gc0_pf"][which]
genome_set_id = [10000, 20000, 60001, 60003][which]

# clean files: IMPORTANT
cmd = """ rm -rf tmp """
run(cmd, shell=True, check=True)


# prepare map files
for chrno in range(1, 15):
    map_fn = f"tmp/map/{chrno}.map"

    Path(map_fn).parent.mkdir(parents=True, exist_ok=True)
    with open(map_fn, "w") as f:
        f.write(f"{chrno} . 0.0 1\n")
        f.write(f"{chrno} . {len_cm} {len_bp}\n")

# prepare sample file
with open("tmp/samples.txt", "w") as f:
    for i in range(0, 1000):
        f.write(f"{i}\n")

# prepare gnome.toml file
genome = {}
len_bp = int(1.0 / r)
# genome file path
# "./r0003_genome.toml"
genome["name"] = f"r{int(r*1e9):04}"
genome["chromsize"] = [
    int(len_bp * 1.2)
] * 14  # make it large to avoid issues due rounding errors
genome["chromnames"] = [f"{i}" for i in range(1, 15)]
genome["gmaps"] = [f"tmp/map/{chrno}.map" for chrno in range(1, 15)]
genome["idx"] = {f"{i}": i - 1 for i in range(1, 15)}

genome_fn = f"./tmp/genome.toml"
with open(genome_fn, "wb") as f:
    tomli_w.dump(genome, f)

root_dir = "."
# copy true ibd
for chrno in range(1, 15):
    # -------------CHANGE HERE AS NECESSARY
    # src_fn = f"res/mp_s00/ibd/tskibd/20000_{chrno}_tskibd.ibd"
    # src_fn = f"res/sp_neu/ibd/tskibd/10000_{chrno}_tskibd.ibd"
    # src_fn = f"res/uk_gc1/ibd/tskibd/60002_{chrno}_tskibd.ibd"
    src_fn = f"{root_dir}/res/{genome_set_name}/ibd/tskibd/{genome_set_id}_{chrno}_tskibd.ibd"
    dst_fn = f"tmp/tskibd/{chrno}.ibd"
    Path(dst_fn).parent.mkdir(parents=True, exist_ok=True)
    # copyfile(src_fn, dst_fn)
    # Path(dst_fn).symlink_to(src_fn)
    Path(dst_fn).hardlink_to(src_fn)

# copy inferred ibd
for p in Path(f"{root_dir}/res/{genome_set_name}/ibd/refinedibd").glob("*"):
    Path(f"tmp/refinedibd/{p.name}").mkdir(exist_ok=True, parents=True)

    for chrno in range(1, 15):
        # -------------CHANGE HERE AS NECESSARY
        # src_fn = f"{p}/20000_{chrno}_hapibd.ibd"
        # src_fn = f"{p}/10000_{chrno}_hapibd.ibd"
        # src_fn = f"{p}/60002_{chrno}_hapibd.ibd"
        src_fn = f"{p}/{genome_set_id}_{chrno}_refinedibd.ibd"
        assert Path(src_fn).exists()
        dst_fn = f"tmp/refinedibd/{p.name}/{chrno}.ibd"
        # copyfile(src_fn, dst_fn)
        # Path(dst_fn).symlink_to(src_fn)
        Path(dst_fn).hardlink_to(src_fn)

# for p in Path("hapibd").glob("*"):
for p in Path("tmp/refinedibd").glob("*"):
    set_id = p.name
    out_prefix = f"tmp/cmp/{set_id}.tsv"
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

    trueibd_dir = "tmp/tskibd"
    infribd_dir = p
    cmd = f"""
    /data/bing/ishare/target/release/ibdutils compare \\
        -g tmp/genome.toml \\
        -s tmp/samples.txt \\
        -S tmp/samples.txt \\
        -f tskibd \\
        -F tskibd \\
        -i {trueibd_dir} \\
        -I {infribd_dir} \\
        -o {out_prefix}
    """
    ret = run(cmd, shell=True, text=True)
    if ret.returncode != 0:
        print("error in running ishare/ibdutils")
        print(f"command: \n {cmd}")
        print(f"\n Error:")
        print(ret.stderr)
        raise ("Error")


# read stats
df_lst = []
keys = ["lod", "len", "scale", "minmac", "window"]
defaults = [4.0, 2.0, 3.16, 20, 40.0]
default_value_str = "lod=4, len=2, scale=3.16, minmac=20, window=40"
for p in Path("tmp/cmp").glob("*.tsv"):
    # lod2.000000_len2.000000_scale2.236068_minmac20_win40.000000
    set_id = p.name.replace(".tsv", "")
    fields = set_id.split("_")
    lod = float(fields[0][3:])
    len = float(fields[1][3:])
    scale = float(fields[2][5:])
    minmac = int(fields[3][6:])
    windows = float(fields[4][3:])
    values = (lod, len, scale, minmac, windows)
    set_str = ""
    for k, v, v2 in zip(keys, defaults, values):
        if k == "lod":
            continue
        if np.abs(v - v2) <= 0.01:
            continue
        if set_str != "":
            set_str += "|"
        set_str += f"{k}{v2:.1f}"
    if set_str == "":
        set_str = "default"
    df = pd.read_csv(p, sep=",")
    df["Set"] = set_str
    df["Lod"] = lod
    df_lst.append(df)
df = pd.concat(df_lst, axis=0)

df.LenBinStart.unique()
lm_lst = df["Lod"].unique()
lm_lst.sort()
set_lst = df["Set"].unique()
set_lst.sort()

cats = ["3", "4", "6", "10", "18", "genome_wide"]


fig = plt.figure(figsize=(12, 6), constrained_layout=True)
gs = fig.add_gridspec(nrows=3, ncols=6, height_ratios=[10, 10, 1])
axes = []
ax_ref = None
for row in [0, 1]:
    ax_row = []
    for col, _ in enumerate(cats):
        if row == 0 and col == 0:
            ax_ref = fig.add_subplot(gs[row, col])
            ax_row.append(ax_ref)
        else:
            ax = fig.add_subplot(gs[row, col])
            ax_row.append(ax)
    axes.append(ax_row)
ax_cbar = fig.add_subplot(gs[2, :])
for col, cat in enumerate(cats):
    # false negative
    # print( df[df.LenBinStart == cat])
    fn = 1 - df[df.LenBinStart == cat].pivot(
        index=["Set"], columns=["Lod"], values=["RateAOverlapByB"]
    )
    # false postive
    fp = 1 - df[df.LenBinStart == cat].pivot(
        index=["Set"], columns=["Lod"], values=["RateBOverlapByA"]
    )
    ax = axes[0][col]
    img = ax.imshow(
        fn, cmap="Blues", interpolation="none", aspect="auto", vmin=0, vmax=1
    )
    for i in range(fn.shape[0]):
        for j in range(fn.shape[1]):
            ax.text(j, i, f"{fn.iloc[i,j]:.2f}", fontsize=7, va="center", ha="center")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticks(np.arange(set_lst.size))
    if col == 0:
        ax.set_yticklabels(set_lst)
    else:
        ax.set_yticklabels([])
    ax = axes[1][col]
    ax.imshow(fp, cmap="Blues", interpolation="none", aspect="auto", vmin=0, vmax=1)
    for i in range(fp.shape[0]):
        for j in range(fp.shape[1]):
            ax.text(j, i, f"{fp.iloc[i,j]:.2f}", fontsize=7, va="center", ha="center")
    ax.set_xticks(np.arange(lm_lst.size))
    ax.set_xticklabels(lm_lst)
    ax.set_yticks(np.arange(set_lst.size))
    if col == 0:
        ax.set_yticklabels(set_lst)
    else:
        ax.set_yticklabels([])
axes[0][0].set_ylabel("\n\n\n\nSet")
axes[1][0].set_ylabel("Set")
plt.colorbar(img, cax=ax_cbar, orientation="horizontal")  # , location="bottom")
for col, _ in enumerate(cats):
    axes[1][col].set_xlabel("Lod")
    axes[0][col].set_title(
        {
            "3": "[3-4) cM",
            "4": "[4-6) cM",
            "6": "[6-10) cM",
            "10": "[10-18) cM",
            "18": "[10-Inf) cM",
            "genome_wide": "genome-wide",
        }[cats[col]],
        fontweight="bold",
    )
axes[0][0].text(
    -1,
    0.5,
    "FN",
    transform=axes[0][0].transAxes,
    in_layout=False,
    fontweight="bold",
    fontsize=10,
)
axes[1][0].text(
    -1,
    0.5,
    "FP",
    transform=axes[1][0].transAxes,
    in_layout=False,
    fontweight="bold",
    fontsize=10,
)
ax_cbar.set_title(default_value_str)

id_desc = f"{demog}"


fig.savefig(f"heatmap_{id_desc}.png", dpi=600)
df.to_csv(f"table_{id_desc}.csv")
