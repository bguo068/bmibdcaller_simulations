from ibdutils.utils.ibdutils import IBD
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import pandas as pd
import tomli_w
from subprocess import run
from shutil import rmtree
import igraph as ig

indir = Path("res/")
indir.mkdir(parents=True, exist_ok=True)

simulations = [
    dict(genome_set_id="40002", label="sp_r0010"),
    dict(genome_set_id="40003", label="sp_r0030"),
    dict(genome_set_id="40004", label="sp_r0100"),
    dict(genome_set_id="40005", label="sp_r0300"),
    dict(genome_set_id="40006", label="sp_r0667"),
    dict(genome_set_id="40007", label="sp_r1000"),
]

ibdcallers = "tskibd hapibd hmmibd".split()

outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)


# ------------------------- IBD comparision -----------------------------------

# prepare map files
len_cm = 100
for sim in simulations:
    label = sim["label"]
    r = int(label[4:]) * 1e-9
    bp_per_cm = 0.01 / r
    len_bp = int(100 * bp_per_cm)

    for chrno in range(1, 15):
        map_fn = f"tmp/{label}/map/{chrno}.map"
        Path(map_fn).parent.mkdir(parents=True, exist_ok=True)
        with open(map_fn, "w") as f:
            f.write(f"{chrno} . 0.0 1\n")
            f.write(f"{chrno} . {len_cm} {len_bp}\n")

    # prepare gnome.toml file
    genome = {}
    # genome file path
    # "./r0003_genome.toml"
    genome["name"] = f"{label}"
    genome["chromsize"] = [
        int(len_bp * 1.001)
    ] * 14  # make it large to avoid issues due rounding errors
    genome["chromnames"] = [f"{i}" for i in range(1, 15)]
    genome["gmaps"] = [f"tmp/{label}/map/{chrno}.map" for chrno in range(1, 15)]
    genome["idx"] = {f"{i}": i - 1 for i in range(1, 15)}

    genome_fn = f"./tmp/{label}/genome.toml"
    with open(genome_fn, "wb") as f:
        tomli_w.dump(genome, f)


# prepare sample file
with open("tmp/samples.txt", "w") as f:
    for i in range(0, 1000):
        f.write(f"{i}\n")


root_dir = "."

for simulation in simulations:
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]

    dir = Path(f"tmp/{label}/tskibd")
    if dir.exists():
        rmtree(dir)
    # link true ibd
    for chrno in range(1, 15):
        # res/mp_s00/ibd/tskibd/20000_11_tskibd.ibd
        src_fn = f"{root_dir}/res/{label}/ibd/tskibd/{genome_set_id}_{chrno}_tskibd.ibd"
        dst_fn = f"tmp/{label}/tskibd/{chrno}.ibd"
        Path(dst_fn).parent.mkdir(parents=True, exist_ok=True)
        Path(dst_fn).hardlink_to(src_fn)

    for caller in ibdcallers:
        if caller == "tskibd":
            caller_ = None
            continue
        if caller == "tpbwt":
            caller_ = "tpbwtibd"
        else:
            caller_ = caller

        dir = Path(f"tmp/{label}/{caller}/")
        if dir.exists():
            rmtree(dir)
        # link inferred ibd
        for chrno in range(1, 15):
            # res/mp_s00/ibd/tskibd/20000_11_tskibd.ibd
            src_fn = f"{root_dir}/res/{label}/ibd/{caller_}/{genome_set_id}_{chrno}_{caller_}.ibd"
            dst_fn = f"tmp/{label}/{caller}/{chrno}.ibd"
            Path(dst_fn).parent.mkdir(parents=True, exist_ok=True)
            Path(dst_fn).hardlink_to(src_fn)

        out_prefix = f"tmp/cmp/{label}/tskibd_vs_{caller}.tsv"
        Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

        # call ishare/ibdutils
        trueibd_dir = f"tmp/{label}/tskibd/"
        infribd_dir = f"tmp/{label}/{caller}/"
        cmd = f"""
        /autofs/chib/toconnor_grp/bing/bmibdcaller_simulations/simulations/r230912/cmp_dwnstrm_ovfilt_nooptimize/ibdutils compare \\
            -g tmp/{label}/genome.toml \\
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


# # read stats
df_lst = []
for simulation in simulations:
    label = simulation["label"]
    for caller in ibdcallers:
        if caller == "tskibd":
            continue
        csv_fn = f"tmp/cmp/{label}/tskibd_vs_{caller}.tsv"
        df = pd.read_csv(csv_fn, sep=",")
        df["Label"] = label
        df["RecomRate"] = int(label[4:]) * 1e-9
        df["Caller"] = caller
        df_lst.append(df)
df = pd.concat(df_lst, axis=0)

df["FN"] = 1 - df["RateAOverlapByB"]
df["FP"] = 1 - df["RateBOverlapByA"]


# %%
nrows = 1
ncols = 2
fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(5, 2),
    constrained_layout=True,
    sharex=True,
    sharey=True,
)

for icaller, caller in enumerate(["hapibd", "hmmibd"]):
    ax = axes[icaller]
    for ilabel, label in enumerate([sim["label"] for sim in simulations]):
        r = int(label[4:]) * 1e-9
        dot_label = f"f={r:.2e}"
        sel = df.LenBinStart == "genome_wide"
        sel2 = (df.Label == label) & (df.Caller == caller)
        df_sel = df[sel & sel2].copy()
        marker = list("osd*hx")[ilabel]
        color = "red blue green purple orange brown cyan".split()[ilabel]
        x, y = df_sel.loc[:, "FN"], df_sel.loc[:, "FP"]
        ax.plot(
            [x],
            [y],
            label=dot_label,
            marker=marker,
            color=color,
            ms=5,
            linestyle="none",
        )
        ax.legend(
            ncol=2, borderpad=0.0, labelspacing=0.1, handletextpad=0.4, handlelength=0.1
        )
        # ax.set_xlim(0, 0.8)
        # ax.set_ylim(0, 0.8)
        ax.set_xlabel("FN")
        if icaller == 0:
            ax.set_ylabel(f"FP")
        ax.set_title(caller)
        # ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        # ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

# %%
# outdir = Path("./analysis_res")
# outdir.mkdir(exist_ok=True, parents=True)
# fig.savefig(f"{outdir}/ibd_cmp_overlap.png", dpi=600)

# ---

totibd_res = {}

# for simulation in simulations:
#     genome_set_id = simulation["genome_set_id"]
#     label = simulation["label"]
#     #
#     for caller in ibdcallers:
#         if caller == "tskibd":
#             continue
#         out_prefix = f"tmp/cmp/{label}/tskibd_vs_{caller}.totibd"
#         df = pd.read_parquet(out_prefix)

#         totibd_res[(label, caller)] = df


# nrows = 4
# ncols = 5
# fig, axes = plt.subplots(
#     nrows=nrows,
#     ncols=ncols,
#     figsize=(10, 8),
#     constrained_layout=True,
#     sharex=True,
#     sharey=True,
# )
# for row, label in enumerate(["sp_neu", "sp_const", "sp_grow", "mp_s00"]):
#     for col, caller in enumerate([c for c in ibdcallers if c != "tskibd"]):
#         df_sel = totibd_res[(label, caller)]
#         df_sel.columns = ["True", "Infer"]
#         df_sel = df_sel.loc[lambda df: ~((df["True"] == 0) & (df.Infer == 0)), :]
#         #
#         ax = axes[row, col]
#         ax.plot(df_sel["True"], df_sel["Infer"], ".", alpha=0.3)
#         ax.plot([2, 1400], [2, 1400])
#         ax.set_yscale("log")
#         ax.set_xscale("log")
#         # ax.legend(
#         #     ncol=2, borderpad=0.0, labelspacing=0.1, handletextpad=0.4, handlelength=0.1
#         # )
#         if row == nrows - 1:
#             ax.set_xlabel("true total IBD")
#         if col == 0:
#             ax.set_ylabel(f"{label}\n\ninferred total IBD")
#         if row == 0:
#             ax.set_title(caller)
#         # ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
#         # ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

# outdir = Path("./analysis_res")
# outdir.mkdir(exist_ok=True, parents=True)
# fig.savefig(f"{outdir}/ibd_cmp_totibd.png", dpi=600)
