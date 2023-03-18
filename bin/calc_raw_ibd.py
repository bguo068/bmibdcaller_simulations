#! /usr/bin/env python3
import argparse
from pathlib import Path
import subprocess
import sys

script_dir = Path(__file__).parent if "__file__" in locals() else Path()
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--chrno", type=int, default=1, help="chromosome number")
parser.add_argument("--bp_per_cm", type=int, default=1000000, help="r = 0.01/bp_per_cm")
parser.add_argument("treeseq_fn", type=str, help="input tree files")
parser.add_argument("--mincm", type=float, help="minimum IBD length in cM")
args = parser.parse_args()
print(args)

# parameters -- simulation
chrno = args.chrno
bp_per_cm = args.bp_per_cm
min_cm = args.mincm
sample_window = int(0.01 * bp_per_cm)
treeseq_fn = args.treeseq_fn

tsk_trueibd_exe = script_dir.parents[1] / "tskibd/build/tskibd"
subprocess.run(
    f"{tsk_trueibd_exe} {chrno} {bp_per_cm} {sample_window} {min_cm} {treeseq_fn}",
    shell=True,
    check=True,
)

# output
"""
{chrno}.log
{chrno}.ibd
{chrno}.map
"""
