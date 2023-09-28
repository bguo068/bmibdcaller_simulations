from ibdutils.utils.ibdutils import IBD
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pathlib import Path
import numpy as np
import pandas as pd
import pickle
import tomli_w
from subprocess import run
from shutil import rmtree
import igraph as ig

indir = Path("./analysis_res/info")
indir.mkdir(parents=True, exist_ok=True)


def pickle_load_file(fn):
    with open(fn, "rb") as f:
        return pickle.load(f)


ibdobj_path_res = pickle_load_file(f"{indir}/ibdobj_path.pkl")
sn_res = pickle_load_file(f"{indir}/signal_noise_ratio.pkl")
np_res = pickle_load_file(f"{indir}/num_peaks.pkl")
ne_res = pickle_load_file(f"{indir}/ne_res.pkl")
ifm_res = pickle_load_file(f"{indir}/infomap_res.pkl")

simulations = [
    dict(genome_set_id="10000", label="sp_neu"),
    dict(genome_set_id="10001", label="sp_s01"),
    dict(genome_set_id="10002", label="sp_s02"),
    dict(genome_set_id="10003", label="sp_s03"),
    dict(genome_set_id="10004", label="sp_g040"),
    dict(genome_set_id="10005", label="sp_g080"),
    dict(genome_set_id="10006", label="sp_g120"),
    dict(genome_set_id="10007", label="sp_o01"),
    dict(genome_set_id="10008", label="sp_o03"),
    dict(genome_set_id="10009", label="sp_o27"),
    dict(genome_set_id="20000", label="mp_s00"),
    dict(genome_set_id="20001", label="mp_s01"),
    dict(genome_set_id="20002", label="mp_s02"),
    dict(genome_set_id="20003", label="mp_s03"),
    dict(genome_set_id="30000", label="sp_rel"),
    dict(genome_set_id="30001", label="mp_rel"),
    dict(genome_set_id="30005", label="sp_const"),
    dict(genome_set_id="30006", label="sp_grow"),
]

ibdcallers = "tskibd hapibd  hmmibd  isorelate  refinedibd  tpbwt".split()

outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)
# ------------------------- plot coverage and peaks --------------------------

for model in ["sp_s03", "mp_s03"]:
    fig, axes = plt.subplots(
        nrows=6,
        ncols=1,
        figsize=(10, 8),
        constrained_layout=True,
        sharex=True,
        sharey=True,
    )
    for icaller, caller in enumerate(ibdcallers):
        ax = axes[icaller]
        if icaller == 0:
            ax.set_title(model)
        IBD.pickle_load(ibdobj_path_res[model, caller]).plot_coverage(
            ax, label=caller, which="xirsfilt", plot_proportions=False
        )
    fig.savefig(f"{outdir}/coverage_{model}.png", dpi=600)


# -------------------------------------------------------------

df = pd.DataFrame(
    [
        {
            "Label": label,
            "Caller": ibdcaller,
            "SNRatio": np.mean(r_lst),
            "Std": np.std(r_lst),
        }
        for (label, ibdcaller), r_lst in sn_res.items()
    ]
)

fig = plt.figure(constrained_layout=True, figsize=(8, 8))
for icaller, label in enumerate(
    ["sp_s01", "mp_s01", "sp_s02", "mp_s02", "sp_s03", "mp_s03"]
):
    df_sel = df[df.Label == label]
    ax = fig.add_subplot(3, 2, icaller + 1)
    ax.bar(df_sel.Caller, df_sel.SNRatio, yerr=df_sel.Std, capsize=5)
    ax.set_ylabel("signal to noise ratio")
    ax.set_title(f"simulation: {label}")
    ax.set_xticklabels(
        ax.get_xticklabels(), rotation=45, rotation_mode="anchor", va="top", ha="right"
    )

fig.savefig(f"{outdir}/signal_noise.png", dpi=600)

# ------------------------------------------------------------------------

df = dict(Label=[], Caller=[], Obs=[], Exp=[])
for (label, caller), (obs, exp) in np_res.items():
    df["Label"].append(label)
    df["Caller"].append(caller)
    df["Obs"].append(obs)
    df["Exp"].append(exp)
df = pd.DataFrame(df)

fig = plt.figure(constrained_layout=True, figsize=(8, 8))
for icaller, label in enumerate(
    ["sp_s01", "mp_s01", "sp_s02", "mp_s02", "sp_s03", "mp_s03"]
):
    # label = "mp_s01"
    df_sel = df[df.Label == label]
    ax = fig.add_subplot(4, 2, icaller + 1)
    ax.bar(df_sel.Caller, df_sel.Obs)
    ax.set_title(f"simulation: {label}")
    ax.set_xticklabels(
        ax.get_xticklabels(), rotation=45, rotation_mode="anchor", va="top", ha="right"
    )
    ax.set_ylabel("no. of peaks observed")
outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)
fig.savefig(f"{outdir}/num_peaks_observed.png", dpi=600)

# --------------------------------------------------------------------------


fig = plt.figure(figsize=(12, 8), constrained_layout=True)
for ilabel, target_label in enumerate(["sp_neu", "sp_const", "sp_grow"]):
    # target_label = "sp_neu"
    target_rm_model = "orig"

    param_df = None
    df_res = {}
    # filter ne dataframes for a given
    for (label, caller, rm_mode), df in ne_res.items():
        if label != target_label:
            continue
        if caller == "param":
            param_df = df
        if rm_mode != target_rm_model:
            continue
        else:
            df_res[caller] = df
            print(label, caller, rm_mode)
    df_res.keys()

    ncallers = len(ibdcallers)
    for icaller, caller in enumerate(ibdcallers):
        ax = fig.add_subplot(3, ncallers, ilabel * ncallers + icaller + 1)
        df_sel = df_res[caller]
        ax.plot(param_df.GEN, param_df.NE, label="Truth", color="k", linestyle="--")
        ax.plot(df_sel.GEN, df_sel.NE, label=caller, color="r")
        ax.fill_between(df_sel.GEN, y1=df_sel.L95, y2=df_sel.U95, fc="r", alpha=0.3)
        ax.set_xlim(0, 100)
        ax.set_ylim(1e2, 1e6)
        ax.set_title(caller)
        ax.legend()
        ax.set_yscale("log")
        if icaller == 0:
            ax.set_ylabel(f"{label.upper()}\n\nNe")
        if ilabel == 2:
            ax.set_xlabel("generations ago")
fig.suptitle("rm_mode = orig")

outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)
fig.savefig(f"{outdir}/ne.png", dpi=600)

# ------------------------- infomap ---------------------------------------


def transform_ifm_df(df):
    """transform infomap df to a 5x5 matrix

    args:
    -----
    df: infomap df (contains Sample, Rank, Population columns)

    columns: community label
    rows: true population label
    cells: number of samples with a given true population label and a given
    community label

    """
    df = (
        df.groupby(["Population", "Rank"])["Sample"]
        .count()
        .unstack()
        .fillna(0)
        .astype(int)
        .iloc[:, :5]
        .copy()
    )
    df = df.sort_values(by=df.index.to_list(), axis=1, ascending=False)
    df.columns = [f"c{i+1}" for i in range(0, 5)]
    df.index = [f"p{i+1}" for i in range(0, 5)]
    return df


def get_adj_rank(df):
    """get adjusted rand index from infomap df"""
    adj_rand = ig.compare_communities(df.Rank, df.Population, method="adjusted_rand")
    return adj_rand


fig = plt.figure(figsize=(12, 8), constrained_layout=True)
fig, axes = plt.subplots(nrows=4, ncols=7, figsize=(12, 8), constrained_layout=True)
for ilabel, label in enumerate(["mp_s00", "mp_s01", "mp_s02", "mp_s03"]):
    adj_rank_lst = []
    for icaller, caller in enumerate(ibdcallers):
        ifm = ifm_res[(label, caller, "orig")]
        adj_rank = get_adj_rank(ifm)
        adj_rank_lst.append(adj_rank)
        ax = axes[ilabel, icaller]
        df = transform_ifm_df(ifm)
        ax.imshow(
            df, cmap="Blues", vmin=0, vmax=200, interpolation="none", aspect="auto"
        )
        ax.set_xticks(np.arange(0, 5))
        ax.set_xticklabels(df.columns)
        ax.set_yticks(np.arange(0, 5))
        ax.set_yticklabels(df.index)
        if ilabel == 0:
            ax.set_title(caller)
        if ilabel == 3:
            ax.set_xlabel("community label")
        if icaller == 0:
            ax.set_ylabel(f"{label}\n\ntrue pop label")
    ax = axes[ilabel, 6]
    ax.bar(np.arange(len(ibdcallers)), adj_rank_lst)
    ax.set_xticks(np.arange(len(ibdcallers)))
    ax.set_xticklabels(ibdcallers, rotation=45, rotation_mode="anchor", ha="right")
    ax.set_ylabel("adjusted rand index")
    ax.set_ylim(0, 1)
    if ilabel == 0:
        ax.set_title("adjusted rand index")
fig.savefig(f"{outdir}/infomap.png", dpi=600)


# ------------------------- IBD comparision -----------------------------------

# prepare map files
bp_per_cm = 15000
len_cm = 100
len_bp = 100 * bp_per_cm
r = 0.01 / bp_per_cm
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
# genome file path
# "./r0003_genome.toml"
genome["name"] = f"r{int(r*1e9):04}"
genome["chromsize"] = [
    int(len_bp * 1.001)
] * 14  # make it large to avoid issues due rounding errors
genome["chromnames"] = [f"{i}" for i in range(1, 15)]
genome["gmaps"] = [f"tmp/map/{chrno}.map" for chrno in range(1, 15)]
genome["idx"] = {f"{i}": i - 1 for i in range(1, 15)}

genome_fn = f"./tmp/genome.toml"
with open(genome_fn, "wb") as f:
    tomli_w.dump(genome, f)


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
for simulation in simulations:
    label = simulation["label"]
    for caller in ibdcallers:
        if caller == "tskibd":
            continue
        csv_fn = f"tmp/cmp/{label}/tskibd_vs_{caller}.tsv"
        df = pd.read_csv(csv_fn, sep=",")
        df["Label"] = label
        df["Caller"] = caller
        df_lst.append(df)
df = pd.concat(df_lst, axis=0)

nrows = 4
ncols = 5
fig, axes = plt.subplots(
    nrows=4,
    ncols=5,
    figsize=(10, 8),
    constrained_layout=True,
    sharex=True,
    sharey=True,
)
for row, label in enumerate(["sp_neu", "sp_const", "sp_grow", "mp_s00"]):
    for col, caller in enumerate([c for c in ibdcallers if c != "tskibd"]):
        df_sel = df[(df.Label == label) & (df.Caller == caller)].copy()
        df_sel["FN"] = 1 - df_sel.RateAOverlapByB
        df_sel["FP"] = 1 - df_sel.RateBOverlapByA
        df_sel = df_sel.set_index("LenBinStart")
        #
        ax = axes[row, col]
        # ax = fig.add_subplot(nrows, ncols, 1 + row * ncols + col)
        bin_map = {
            "3": "[3-4)",
            "4": "[4-6)",
            "6": "[6-10)",
            "10": "[10-18)",
            "18": "[18-Inf)",
            "genome_wide": "gw",
        }
        for i, b in enumerate(["3", "4", "6", "10", "18", "genome_wide"]):
            marker = list("osd*hx")[i]
            color = "red blue green purple orange brown cyan".split()[i]
            x, y = df_sel.loc[b, "FN"], df_sel.loc[b, "FP"]
            ax.plot(
                [x],
                [y],
                label=bin_map[b],
                marker=marker,
                color=color,
                ms=5,
                linestyle="none",
            )
        ax.legend(
            ncol=2, borderpad=0.0, labelspacing=0.1, handletextpad=0.4, handlelength=0.1
        )
        ax.set_xlim(0, 0.8)
        ax.set_ylim(0, 0.8)
        if row == nrows - 1:
            ax.set_xlabel("FN")
        if col == 0:
            ax.set_ylabel(f"{label}\n\nFP")
        if row == 0:
            ax.set_title(caller)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)
fig.savefig(f"{outdir}/ibd_cmp_overlap.png", dpi=600)

# ---

totibd_res = {}

for simulation in simulations:
    genome_set_id = simulation["genome_set_id"]
    label = simulation["label"]
    #
    for caller in ibdcallers:
        if caller == "tskibd":
            continue
        out_prefix = f"tmp/cmp/{label}/tskibd_vs_{caller}.totibd"
        df = pd.read_parquet(out_prefix)

        totibd_res[(label, caller)] = df


nrows = 4
ncols = 5
fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(10, 8),
    constrained_layout=True,
    sharex=True,
    sharey=True,
)
for row, label in enumerate(["sp_neu", "sp_const", "sp_grow", "mp_s00"]):
    for col, caller in enumerate([c for c in ibdcallers if c != "tskibd"]):
        df_sel = totibd_res[(label, caller)]
        df_sel.columns = ["True", "Infer"]
        df_sel = df_sel.loc[lambda df: ~((df["True"] == 0) & (df.Infer == 0)), :]
        #
        ax = axes[row, col]
        ax.plot(df_sel["True"], df_sel["Infer"], ".", alpha=0.3)
        ax.plot([2, 1400], [2, 1400])
        ax.set_yscale("log")
        ax.set_xscale("log")
        # ax.legend(
        #     ncol=2, borderpad=0.0, labelspacing=0.1, handletextpad=0.4, handlelength=0.1
        # )
        if row == nrows - 1:
            ax.set_xlabel("true total IBD")
        if col == 0:
            ax.set_ylabel(f"{label}\n\ninferred total IBD")
        if row == 0:
            ax.set_title(caller)
        # ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        # ax.yaxis.set_major_locator(MaxNLocator(nbins=4))

outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)
fig.savefig(f"{outdir}/ibd_cmp_totibd.png", dpi=600)


# -------------------------- binned population-level total IBD ----------------
def get_bin_pop_total_ibd(label, caller):
    ibd = IBD.pickle_load(ibdobj_path_res[(label, caller)])
    cm = ibd.calc_ibd_length_in_cm()
    bins = np.arange(2, 100, 0.05)
    bin_center = (bins[:-1] + bins[1:]) / 2
    bin_count, _ = np.histogram(cm, bins=bins)
    return (bin_center, bin_center * bin_count)


nrows = 4
ncols = 5
fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(10, 8),
    constrained_layout=True,
    sharex=True,
    sharey=True,
)
for row, label in enumerate(["sp_neu", "sp_const", "sp_grow", "mp_s00"]):
    bin_center, x0 = get_bin_pop_total_ibd(label, "tskibd")
    for col, caller in enumerate([c for c in ibdcallers if c != "tskibd"]):
        # true
        x = x0.copy()
        # infer
        _, y = get_bin_pop_total_ibd(label, caller)
        sel1 = (x == 0) & (y == 0)
        sel2 = (x == 0) & (y != 0)
        sel3 = (x != 0) & (y == 0)
        sel4 = (x != 0) & (y != 0)
        x[sel2] = 1
        y[sel3] = 1
        ax = axes[row, col]
        ax.plot(bin_center[sel4], y[sel4] / x[sel4], ".", color="blue", alpha=0.3)
        ax.plot(
            bin_center[sel2 | sel3],
            y[sel2 | sel3] / x[sel2 | sel3],
            ".",
            color="red",
            alpha=0.3,
        )
        ax.axhline(y=1.0, color="r", linestyle="--")
        ax.set_yscale("log")
        ax.set_xscale("log")
        if row == nrows - 1:
            ax.set_xlabel("true total IBD")
        if col == 0:
            ax.set_ylabel(f"{label}\n\ninferred total IBD")
        if row == 0:
            ax.set_title(caller)
        # ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        # ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
fig.suptitle("Population-level total length of IBD of 0.05cM length windows")

outdir = Path("./analysis_res")
outdir.mkdir(exist_ok=True, parents=True)
fig.savefig(f"{outdir}/ibd_cmp_pop_level_totibd.png", dpi=600)
