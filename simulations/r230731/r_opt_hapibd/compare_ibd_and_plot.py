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
p.add_argument("id", type=int, choices=list(range(12)))
which = p.parse_args().id

# clean files
cmd = """
rm -rf cmp hapibd map tskibd samples.txt genome.toml
"""
run(cmd, shell=True, check=True)

# sample file (shared by all comparision as nsam is the same)
with open("samples.txt", "w") as f:
    for i in range(0, 1000):
        f.write(f"{i}\n")


# prepare genome files
# run ishare
# analyze ishare results

# -------------CHANGE HERE AS NECESSARY
r = 0.01 / ([15000] * 4 + [1000000] * 4 + [15000] * 4)[which]
len_bp = ([15000 * 100] * 4 + [1000000 * 60] * 4 + [15000 * 100] * 4)[which]
len_cm = ([100.0] * 4 + [60.0] * 4 + [100.0] * 4)[which]


# prepare map files
for chrno in range(1, 15):
    map_fn = f"map/{chrno}.map"

    Path(map_fn).parent.mkdir(parents=True, exist_ok=True)
    with open(map_fn, "w") as f:
        f.write(f"{chrno} . 0.0 1\n")
        f.write(f"{chrno} . {len_cm} {len_bp}\n")

# ---------------------------------------
# gnome
genome = {}
len_bp = int(1.0 / r)
# genome file path
# "./r0003_genome.toml"
genome["name"] = f"r{int(r*1e9):04}"
genome["chromsize"] = [
    int(len_bp * 1.2)
] * 14  # make it large to avoid issues due rounding errors
genome["chromnames"] = [f"{i}" for i in range(1, 15)]
genome["gmaps"] = [f"map/{chrno}.map" for chrno in range(1, 15)]
genome["idx"] = {f"{i}": i - 1 for i in range(1, 15)}

genome_fn = f"./genome.toml"
with open(genome_fn, "wb") as f:
    tomli_w.dump(genome, f)

min_mac = [20, 200, 20, 200, 20, 200, 20, 200, 20, 200, 20, 200][which]
genome_set_name = [
    "sp_neu",
    "sp_neu",
    "mp_s00",
    "mp_s00",
    "uk_gc0",
    "uk_gc0",
    "uk_gc1",
    "uk_gc1",
    "uk_gc0_pf",
    "uk_gc0_pf",
    "uk_gc1_pf",
    "uk_gc1_pf",
][which]
genome_set_id = [
    10000,
    10000,
    20000,
    20000,
    60001,
    60001,
    60002,
    60002,
    60003,
    60003,
    60004,
    60004,
][which]
# copy true ibd
for chrno in range(1, 15):
    # -------------CHANGE HERE AS NECESSARY
    # src_fn = f"res/mp_s00/ibd/tskibd/20000_{chrno}_tskibd.ibd"
    # src_fn = f"res/sp_neu/ibd/tskibd/10000_{chrno}_tskibd.ibd"
    # src_fn = f"res/uk_gc1/ibd/tskibd/60002_{chrno}_tskibd.ibd"
    src_fn = f"res/{genome_set_name}/ibd/tskibd/{genome_set_id}_{chrno}_tskibd.ibd"
    dst_fn = f"tskibd/{chrno}.ibd"
    Path(dst_fn).parent.mkdir(parents=True, exist_ok=True)
    copyfile(src_fn, dst_fn)

# copy inferred ibd
# -------------CHANGE HERE AS NECESSARY
# for p in Path("res/mp_s00/ibd/hapibd").glob("*_mac200"):
for p in Path(f"res/{genome_set_name}/ibd/hapibd").glob(f"*_mac{min_mac:03d}"):
    Path(f"hapibd/{p.name}").mkdir(exist_ok=True, parents=True)

    for chrno in range(1, 15):
        # -------------CHANGE HERE AS NECESSARY
        # src_fn = f"{p}/20000_{chrno}_hapibd.ibd"
        # src_fn = f"{p}/10000_{chrno}_hapibd.ibd"
        # src_fn = f"{p}/60002_{chrno}_hapibd.ibd"
        src_fn = f"{p}/{genome_set_id}_{chrno}_hapibd.ibd"
        dst_fn = f"hapibd/{p.name}/{chrno}.ibd"
        copyfile(src_fn, dst_fn)

# for p in Path("hapibd").glob("*"):
for p in Path("hapibd").glob(f"*_mac{min_mac:03d}"):
    set_id = p.name
    out_prefix = f"cmp/{set_id}.tsv"
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

    trueibd_dir = "tskibd"
    infribd_dir = p
    cmd = f"""
    /data/bing/ishare/target/release/ibdutils compare \\
        -g genome.toml \\
        -s samples.txt \\
        -S samples.txt \\
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
for p in Path("cmp").glob("*.tsv"):
    set_id = p.name.replace(".tsv", "")
    mg = int(set_id.split("_")[3][2:])
    mm = int(set_id.split("_")[4][2:])
    df = pd.read_csv(p, sep=",")
    df["MaxGap"] = mg
    df["MinMarkers"] = mm
    df_lst.append(df)
df = pd.concat(df_lst, axis=0)

df.LenBinStart.unique()
markers = df["MinMarkers"].unique()
markers.sort()
maxgaps = df["MaxGap"].unique()
maxgaps.sort()

cats = ["3", "4", "6", "10", "18", "genome_wide"]


fig = plt.figure(figsize=(12, 5), constrained_layout=True)
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
    fn = 1 - df[df.LenBinStart == cat].pivot(
        index=["MaxGap"], columns=["MinMarkers"], values=["RateAOverlapByB"]
    )
    # false postive
    fp = 1 - df[df.LenBinStart == cat].pivot(
        index=["MaxGap"], columns=["MinMarkers"], values=["RateBOverlapByA"]
    )
    ax = axes[0][col]
    img = ax.imshow(fn, cmap="Blues", interpolation="none", aspect=True, vmin=0, vmax=1)
    for i in range(fn.shape[0]):
        for j in range(fn.shape[1]):
            ax.text(j, i, f"{fn.iloc[i,j]:.2f}", fontsize=7, va="center", ha="center")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticks(np.arange(maxgaps.size))
    if col == 0:
        ax.set_yticklabels(maxgaps)
    else:
        ax.set_yticklabels([])
    ax = axes[1][col]
    ax.imshow(fp, cmap="Blues", interpolation="none", aspect=True, vmin=0, vmax=1)
    for i in range(fp.shape[0]):
        for j in range(fp.shape[1]):
            ax.text(j, i, f"{fp.iloc[i,j]:.2f}", fontsize=7, va="center", ha="center")
    ax.set_xticks(np.arange(markers.size))
    ax.set_xticklabels(markers)
    ax.set_yticks(np.arange(maxgaps.size))
    if col == 0:
        ax.set_yticklabels(maxgaps)
    else:
        ax.set_yticklabels([])
axes[0][0].set_ylabel("\n\n\n\nMaxGap")
axes[1][0].set_ylabel("MaxGap")
plt.colorbar(img, cax=ax_cbar, orientation="horizontal")  # , location="bottom")
for col, _ in enumerate(cats):
    axes[1][col].set_xlabel("MinMarkers")
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
    -0.7,
    0.5,
    "FN",
    transform=axes[0][0].transAxes,
    in_layout=False,
    fontweight="bold",
    fontsize=10,
)
axes[1][0].text(
    -0.7,
    0.5,
    "FP",
    transform=axes[1][0].transAxes,
    in_layout=False,
    fontweight="bold",
    fontsize=10,
)

demog = [
    "sp",
    "sp",
    "mp",
    "mp",
    "uknogc",
    "uknogc",
    "ukgc",
    "ukgc",
    "uknogcpf",
    "uknogcpf",
    "ukgcpf",
    "ukgcpf",
][which]
id_desc = f"minmarker_maxgap_{demog}_minmac{min_mac}"


fig.savefig(f"heatmap_{id_desc}.png", dpi=600)
df.to_csv(f"table_{id_desc}.csv")
