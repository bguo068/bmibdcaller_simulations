#! /usr/bin/env python3
from pathlib import Path
from shutil import copyfile, rmtree
from subprocess import run

import matplotlib.pyplot as plt
import pandas as pd
import tomli_w


# sample file (shared by all comparision as nsam is the same)
with open("samples.txt", "w") as f:
    for i in range(0, 1000):
        f.write(f"{i}\n")


# prepare genome files
# run ishare
# analyze ishare results

# %%
for x_ in ["hapibd", "hmmibd", "refinedibd"]:
    ibdcallers = ["tskibd", x_]
    n_valid_chr_dict = {}
    model = ["single", "multiple"][1]
    print(ibdcallers, model)
    for ir, r in enumerate([3, 10, 30, 100, 300, 667, 1000]):
        if Path(f"r{r:04}/map").exists():
            rmtree(f"r{r:04d}/map")
        if Path(f"r{r:04}/ibd").exists():
            rmtree(f"r{r:04d}/ibd")
        genome = {}
        len_bp = int(1.0 / r * 1e9)

        # prepare map files
        for chrno in range(1, 15):
            map_fn = f"r{r:04}/map/{chrno}.map"

            Path(map_fn).parent.mkdir(parents=True, exist_ok=True)
            with open(map_fn, "w") as f:
                f.write(f"{chrno} . 0.0 1\n")
                f.write(f"{chrno} . 100.0 {len_bp}\n")

        # genome file path
        # "./r0003_genome.toml"
        genome["name"] = f"r{r:04}"
        genome["chromsize"] = [
            int(len_bp * 1.2)
        ] * 14  # make it large to avoid issues due rounding errors
        genome["chromnames"] = [f"{i}" for i in range(1, 15)]
        genome["gmaps"] = [f"r{r:04}/map/{chrno}.map" for chrno in range(1, 15)]
        genome["idx"] = {f"{i}": i - 1 for i in range(1, 15)}

        genome_fn = f"./r{r:04d}_genome.toml"
        with open(genome_fn, "wb") as f:
            tomli_w.dump(genome, f)

        # copy/reorganize ibd files
        n_valid_chr = 0
        for chrno in range(1, 15):
            n_exist = 0
            modelabbr = "s" if model == "single" else "m"
            setidabbr = "4" if model == "single" else "5"
            for ibdcaller in ibdcallers:
                src_fn = f"../res/{modelabbr}p_r{r:04}/ibd/{ibdcaller}/{setidabbr}000{ir+1}_{chrno}_{ibdcaller}.ibd"
                if Path(src_fn).exists():
                    n_exist += 1
            if n_exist != 2:
                continue
            for ibdcaller in ibdcallers:
                src_fn = f"../res/{modelabbr}p_r{r:04}/ibd/{ibdcaller}/{setidabbr}000{ir+1}_{chrno}_{ibdcaller}.ibd"
                dst_fn = f"r{r:04}/ibd/{ibdcaller}/{chrno}.ibd"
                Path(dst_fn).parent.mkdir(parents=True, exist_ok=True)
                copyfile(src_fn, dst_fn)
            n_valid_chr += 1
        n_valid_chr_dict[r] = n_valid_chr

    for r in [10, 30, 100, 300, 667, 1000]:
        if n_valid_chr_dict[r] < 10:
            continue

        genome_fn = f"./r{r:04d}_genome.toml"

        ibd_dirs = []
        for ibdcaller in ibdcallers:
            ibd_dirs.append(f"r{r:04}/ibd/{ibdcaller}")

        out_prefix = f"r{r:04}/cmp_{ibdcallers[0]}_{ibdcallers[1]}_{model}.tsv"
        cmd = f"""
        /data/bing/ishare/target/release/ibdutils compare \\
            -g {genome_fn} \\
            -s samples.txt  \\
            -S samples.txt  \\
            -f tskibd \\
            -F tskibd \\
            -i {ibd_dirs[0]}  \\
            -I {ibd_dirs[1]} \\
            -o {out_prefix}
            """
        print(r)
        ret = run(cmd, shell=True, text=True)
        if ret.returncode != 0:
            print("error in running ishare/ibdutils")
            print(f"command: \n {cmd}")
            print("\n Error:")
            print(ret.stderr)
            raise ("Error")

    # read stats
    df_lst = []
    for r in [10, 30, 100, 300, 667, 1000]:
        # TODO: skip 3 for now became the simulation for r==3 is done yet
        if r not in n_valid_chr_dict or n_valid_chr_dict[r] < 10:
            continue

        out_prefix = f"r{r:04}/cmp_{ibdcallers[0]}_{ibdcallers[1]}_{model}.tsv"
        df = pd.read_csv(out_prefix, sep=",")
        df["Rec"] = r
        df_lst.append(df)
    df = pd.concat(df_lst, axis=0)

    fig, axes = plt.subplots(nrows=3, figsize=(6, 9), constrained_layout=True)
    df3 = df[df.LenBinStart == "genome_wide"].sort_values("Rec")
    for r in [3, 10, 30, 100, 300, 667, 1000]:
        if r not in n_valid_chr_dict or n_valid_chr_dict[r] < 10:
            continue
        df2 = df[(~(df.LenBinStart == "genome_wide")) & (df.Rec == r)].copy()
        df2["LenBinStart"] = df2["LenBinStart"].astype(int)
        ax = axes[0]
        ax.plot(
            df2["LenBinStart"],
            1 - df2["RateAOverlapByB"],
            marker=".",
            ms=10,
            label=f"r={r/1e9:.1e}",
        )
        ax.set_xlabel("Length Bin")
        ax.set_ylabel("False neg")
        ax.legend()
        ax = axes[1]
        ax.plot(
            df2["LenBinStart"],
            1 - df2["RateBOverlapByA"],
            marker=".",
            ms=10,
            label=f"r={r/1e9:.1e}",
        )
        ax.set_xlabel("Lengh Bin")
        ax.set_ylabel("False pos")
        ax.legend()
    ax = axes[2]
    for lenstart, a_by_b, b_by_a, r in df3.itertuples(index=None):
        ax.plot(
            [1 - a_by_b],
            [1 - b_by_a],
            marker=".",
            linestyle="none",
            ms=10,
            label=f"r={r/1e9:.2e}",
        )
    ax.set_ylabel("False posive")
    ax.set_xlabel("False negative")
    ax.set_title("pair-wise overlap")
    ax.legend()
    axes[0].set_ylim(0, 1)
    axes[1].set_ylim(0, 1)
    axes[2].set_xlim(0, 1)
    axes[2].set_ylim(0, 1)
    fig.savefig(f"cmp_{ibdcallers[0]}_{ibdcallers[1]}_{model}.png", dpi=600)

# %%
