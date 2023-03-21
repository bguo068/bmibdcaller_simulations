#! /usr/bin/env python3

from argparse import ArgumentParser

from ibdutils.utils.ibdutils import IBD, Genome, IbdComparator

p = ArgumentParser()
p.add_argument("--true_ibd_lst", nargs="+", type=str, required=True)
p.add_argument("--inferred_ibd_lst", nargs="+", type=str, required=True)
p.add_argument("--r", type=float, required=True)
p.add_argument("--seqlen", type=int, required=True)
p.add_argument("--nchroms", type=int, required=True)
p.add_argument("--ibdcaller", type=str, required=True, help="caller for infered ibd")
p.add_argument("--genome_set_id", type=int, required=True)

args = p.parse_args()
nchroms = args.nchroms
r = args.r
seqlen = args.seqlen

assert len(args.true_ibd_lst) == nchroms
assert len(args.inferred_ibd_lst) == nchroms

genome = Genome.get_genome_simple_simu(r, nchroms, seqlen)

ibd1 = IBD(genome, label="tskibd")
ibd1.read_ibd(args.true_ibd_lst)

ibd2 = IBD(genome, label=args.ibdcaller)
ibd2.read_ibd(args.inferred_ibd_lst)

"""
ibd dataframe format
Id1  Id2    Start      End  Ancestor  Tmrca  HasMutation  Chromosome
"""


cmp = IbdComparator(ibd1, ibd2)

cmp.compare()

ofn = f"{args.genome_set_id}_{args.ibdcaller}.ibdcmpobj.gz"

cmp.pickle_dump(ofn)
