#! /usr/bin/env python3
import argparse
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parents[2] / "src"))
# print(__file__)
# print(sys.path)
import smsimu

script_dir = Path(__file__).parent if "__file__" in locals() else Path()
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--chrno", type=int, default=1, help="chromosome number")
parser.add_argument("--s", type=float, default=0.1, help="selection coefficient")
parser.add_argument("--h", type=float, default=0.5, help="dominance coefficient")
parser.add_argument("--selpos_0_1", type=float, default=0.33, help="selpos: 0~1")
parser.add_argument("--num_origins", type=int, default=1, help="dac under sel")
parser.add_argument("--s_start_g", type=int, default=50, help="selStartG")
parser.add_argument("--u", type=float, default=1e-8, help="only sim anc, u ignored")
parser.add_argument("--bp_per_cm", type=int, default=1000000, help="r = 0.01/bp_per_cm")
parser.add_argument("--seqlen_in_cm", type=int, default=100, help="chrlen in cM")
parser.add_argument(
    "--ne_change_start_g", type=int, default=200, help="Ne change time in gen ago)"
)
parser.add_argument("--N", type=int, default=10000, help="ancient Ne")
parser.add_argument("--N0", type=int, default=10000, help="Ne at present")
parser.add_argument("--nsam", type=int, default=1000, help="no. hap sampled")
parser.add_argument(
    "--test",
    action="store_true",
    dest="test",
    help="when this is set, nsam values is replaced with 50 and seqlen with 100cM for quick testing",
)
parser.add_argument(
    "--sim_related",
    type=int,
    choices=[0, 1],
    help="0: normal simulation; 1: simulate highly related",
)
args = parser.parse_args()
print(args)

# parameters -- simulation
sim_related = args.sim_related
chrno = args.chrno
nsam = args.nsam
s_start_g = args.s_start_g
num_origins = args.num_origins
h = args.h
bp_per_cm = args.bp_per_cm
r = 0.01 / bp_per_cm
u = args.u
N = args.N
N0 = args.N0
ne_change_start_g = args.ne_change_start_g
s = args.s
seqlen = args.seqlen_in_cm * bp_per_cm

# parameters -- find true ibd
min_cm = 2.0
sample_window = int(0.01 * bp_per_cm)
min_bp = min_cm * bp_per_cm
remove_hbd = True
min_tmrca = 1.5
out_ibd = f"{chrno}.ibd"

# if test using a small number for nsam and seqlen
if args.test:
    nsam = 50
    seqlen = 20 * bp_per_cm

selpos_bp = seqlen * args.selpos_0_1

# setting parameters for simulator wrapper
simulator = smsimu.SlimMsprimeSimulator(
    seqlen=seqlen,
    r=r,
    u=u,
    N=N,
    N0=N0,
    ne_change_g=ne_change_start_g,
    nsam=nsam,
    nrep=1,
    h=h,
    s=s,
    selpos_bp=selpos_bp,
    s_start_g=s_start_g,
    num_origins=num_origins,
    sim_related=sim_related,
)

# first do slim simulation and tree sequence recording, then using
# msprime/pyslim to recapitate (finish coalescence)
ts, slim_restart_count = simulator.simulate_a_chromosome(
    idx=chrno,
    slim_seed=chrno,
    recapitate_seed=chrno * chrno,
    skip_sim_mutation=True,  # allow to only analyze mutation under selection from slim
)

Path(f"{chrno}.restart_count").write_text(f"{slim_restart_count}")
simulator.true_ne_df.to_csv(f"{chrno}.true_ne", sep="\t", index=None)
simulator.daf_df.to_csv(f"{chrno}.daf", index=None, sep="\t")

treeseq_fn = f"{chrno}.trees"
ts.dump(treeseq_fn)

# output files:
"""
{chrno}.trees
tmp_slim_out_{chrno}.trees
{chrno}.daf
{chrno}.true_ne
{chrno}.restart_count
"""
