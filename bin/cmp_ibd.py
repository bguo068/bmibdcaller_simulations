#! /usr/bin/env python3
import tomli_w
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from subprocess import run
from shutil import copyfile, rmtree
import numpy as np
import argparse

"""
# move ibdutils to bin folder
# install tomli_w

num_chrs // params.nchroms
inferred_ibd_dir
true_ibd_dir
args.seqlen_cm
args.r

# steps
- preprare files: a toml file,  map files and a sample list file

# run ibdutils tools

# output overlapping rate table

# output total ibd parquet file
"""

p = argparse.ArgumentParser()
p.add_argument("--id", type=str, required=True)
p.add_argument("--r", type=float, required=True)
p.add_argument("--nchroms", type=int, required=True)
p.add_argument("--seqlen_bp", type=int, required=True)
p.add_argument("--trueibd_dir", type=str, required=True)
p.add_argument("--inferredibd_dir", type=str, required=True)
args = p.parse_args()

# groups
r = args.r
len_bp = args.seqlen_bp
len_cm = 100 * r * len_bp
set_id = args.id
trueibd_dir = args.trueibd_dir
inferredibd_dir = args.inferredibd_dir

# prepare map files
for chrno in range(1, 15):
    map_fn = f"map/{chrno}.map"

    Path(map_fn).parent.mkdir(parents=True, exist_ok=True)
    with open(map_fn, "w") as f:
        f.write(f"{chrno} . 0.0 1\n")
        f.write(f"{chrno} . {len_cm} {len_bp}\n")

# prepare sample file
with open("samples.txt", "w") as f:
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
genome["gmaps"] = [f"map/{chrno}.map" for chrno in range(1, 15)]
genome["idx"] = {f"{i}": i - 1 for i in range(1, 15)}

genome_fn = f"genome.toml"
with open(genome_fn, "wb") as f:
    tomli_w.dump(genome, f)


out_prefix = f"cmp.ovcsv"
Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

# overlap, totibd analysis
cmd = f"""
ibdutils --version
ibdutils compare \\
    -g genome.toml \\
    -s samples.txt \\
    -S samples.txt \\
    -f tskibd \\
    -F tskibd \\
    -i {trueibd_dir} \\
    -I {inferredibd_dir} \\
    -o {out_prefix}
"""
ret = run(cmd, shell=True, text=True)

print(cmd)
print(set_id)

if ret.returncode != 0:
    print("error in running ishare/ibdutils")
    print(f"\n Error:")
    print(ret.stderr)
    raise ("Error")
