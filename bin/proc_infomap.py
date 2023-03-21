#! /usr/bin/env python3

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import ibdutils.utils.ibdutils as ibdutils

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs=14, required=True)
pa.add_argument("--genome_set_id", type=int, required=True)
pa.add_argument("--r", type=float, default=0.01 / 15000)
pa.add_argument("--seqlen", type=int, default=100 * 15000)
pa.add_argument("--nchroms", type=int, default=14)
args = pa.parse_args()


# read ibd

genome = ibdutils.Genome.get_genome_simple_simu(args.r, args.nchroms, args.seqlen)
ibd = ibdutils.IBD(genome=genome, label="orig")
# remove sample name suffix (haploid genome ratio)
# ibd.read_ibd(ibd_fn_lst=args.ibd_files, rm_sample_name_suffix=True)
ibd.read_ibd(ibd_fn_lst=args.ibd_files)


# output files:
of_ifm_orig_ibdobj = f"{args.genome_set_id}_ifm_orig.ibdobj.gz"
of_ifm_rmpeaks_ibdobj = f"{args.genome_set_id}_ifm_rmpeaks.ibdobj.gz"


# remove highly relatedness samples
mat = ibd.make_ibd_matrix()
unrelated_samples = ibd.get_unrelated_samples(ibdmat=mat)
ibd.subset_ibd_by_samples(subset_samples=unrelated_samples)


# calculate coverage and remove peaks
ibd.calc_ibd_cov()
ibd.find_peaks()

ibd2 = ibd.duplicate("rmpeaks")
ibd2.remove_peaks()

ibd.pickle_dump(of_ifm_orig_ibdobj)
ibd2.pickle_dump(of_ifm_rmpeaks_ibdobj)
