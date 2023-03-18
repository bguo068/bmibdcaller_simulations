#! /usr/bin/env python3
import pyslim
import msprime
import tskit
import subprocess
from pathlib import Path
import numba
import pandas as pd
import numpy as np
import sys
import argparse

sys.path.append(str(Path(__file__).parents[1] / "src"))
print(sys.path)

import trueibd

parser = argparse.ArgumentParser()
parser.add_argument("--chrno", type=int, required=True)
parser.add_argument("--bp_per_cm", type=int, default=1000000)
parser.add_argument("--test", action="store_true", dest="test")
args = parser.parse_args()

# parameters
chrno = args.chrno
bp_per_cm = args.bp_per_cm
N = 10000
r = 0.01 / bp_per_cm
u = 1e-8
seqlen = 100 * bp_per_cm  ## TEST 20 * bp_per_cm
nsam = 1000  ## TEST 50
if args.test:
    seqlen = 20 * bp_per_cm  ## TEST 20 * bp_per_cm
    nsam = 50  ## TEST 50

min_cm = 2
sample_window = int(0.01 * bp_per_cm)
remove_hbd = True
merge_ibd = True
min_tmrca = 1.5
min_bp = bp_per_cm * min_cm
slim_total_generations = 200
slim_tree = Path(f"tmp_{chrno}.trees")
slim_script = Path(f"tmp_{chrno}.slim")
out_tree = Path(f"{chrno}.trees")
out_ibd = Path(f"{chrno}.ibd")

# slim script
slim_script_str = f"""
// adapted_frrecipe_17.1.trees
initialize() {{
	initializeTreeSeq();
	initializeMutationRate(0);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, {seqlen}-1);
	initializeRecombinationRate({r});
}}
1 {{
	sim.addSubpop("p1", {N});
}}
{slim_total_generations} late() {{
	sim.treeSeqOutput("{slim_tree}");
}}
"""
slim_script.write_text(slim_script_str)

# run slim simulation
subprocess.run(["slim", "-s", f"{chrno}", str(slim_script)], check=True)
#
# load trees
orig_ts = pyslim.load(str(slim_tree))

# recapitate
rts = pyslim.recapitate(
    orig_ts, recombination_rate=r, ancestral_Ne=N, random_seed=chrno
)

# simplification
np.random.seed(chrno)
alive_inds = rts.individuals_alive_at(0)
keep_indivs = np.random.choice(alive_inds, nsam // 2, replace=False)
keep_nodes = []
for i in keep_indivs:
    keep_nodes.extend(rts.individual(i).nodes)

ts = rts.simplify(keep_nodes, keep_input_roots=False)

# find true ibd
df_ibd_raw = trueibd.find_true_ibd_from_ts(ts, sample_window, min_bp=2.0 * bp_per_cm)
print(df_ibd_raw.shape)
df_ibd, df_ibd_unmerged = trueibd.make_diploid_filterred_merged_ibd(
    df_ibd_raw, bp_per_cm, chrno=chrno, return_unmerge_ibd=True
)

# output ibd
df_ibd_raw.to_csv(out_ibd.with_suffix(".ibd_raw"), sep="\t", index=None)
df_ibd_unmerged.to_csv(out_ibd.with_suffix(".ibd_unmerged"), sep="\t", index=None)
df_ibd.to_csv(out_ibd, sep="\t", index=None, header=None)

df_map = trueibd.make_plink_map(chrno, bp_per_cm, seqlen)
df_map.to_csv(f"{chrno}.map")
