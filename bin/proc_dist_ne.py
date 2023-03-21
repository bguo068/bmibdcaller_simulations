#! /usr/bin/env python3

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path

import ibdutils.runner.ibdne as ibdne
import ibdutils.utils.ibdutils as ibdutils

ibd_jar_default = str(Path(__file__).parents[1] / "lib/ibdne.23Apr20.ae9.jar")

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs="+", required=True)
pa.add_argument("--genome_set_id", type=int, required=True)
pa.add_argument("--ibdne_jar", type=str, default=ibd_jar_default)
pa.add_argument("--ibdne_mincm", type=float, default=4)
pa.add_argument("--ibdne_minregion", type=float, default=10)
pa.add_argument("--r", type=float, default=0.01 / 15000)
pa.add_argument("--seqlen", type=int, default=100 * 15000)
pa.add_argument("--nchroms", type=int, default=14)

args = pa.parse_args()

assert args.ibdne_jar is not None

genome_set_id = args.genome_set_id
ibdne_mincm = args.ibdne_mincm
ibdne_minregion = args.ibdne_minregion

# create genome obj
genome = ibdutils.Genome.get_genome_simple_simu(args.r, args.nchroms, args.seqlen)

# create ibd obj
ibd = ibdutils.IBD(genome=genome, label=f"{genome_set_id}_orig")

# read ibd from files each containing IBD segment info for a chromosome
# and remove sample name suffix (such as "~0.9-0.2") that was used
# for carrying haploid genome ratios
# ibd.read_ibd(ibd_fn_lst=args.ibd_files, rm_sample_name_suffix=True)
ibd.read_ibd(ibd_fn_lst=args.ibd_files)

# ---- save origin IBD obj for distribution analyses
of_ibddist_ibdobj = f"{genome_set_id}_ibddist.ibdobj.gz"
ibd.pickle_dump(of_ibddist_ibdobj)


# ---- save IBD with coverage for all samples (haploid)
ibd_all = ibd.duplicate()
ibd_all.filter_ibd_by_length(min_seg_cm=ibdne_mincm)
ibd_all.calc_ibd_cov()
ibd_all.find_peaks()
ibd_all.pickle_dump(f"{genome_set_id}_orig_all.ibdcov.ibdobj.gz")

# remove highly related samples
mat = ibd.make_ibd_matrix()
unrelated_samples = ibd.get_unrelated_samples(ibdmat=mat)
ibd.subset_ibd_by_samples(subset_samples=unrelated_samples)

# ---- save ibd coverage (haploid) for unrelated samples
ibd_unrel = ibd.duplicate()
ibd_unrel.filter_ibd_by_length(min_seg_cm=ibdne_mincm)
ibd_unrel.calc_ibd_cov()
ibd_unrel.find_peaks()
ibd_unrel.pickle_dump(f"{genome_set_id}_orig_unrel.ibdcov.ibdobj.gz")

# remove short segments
ibd.filter_ibd_by_length(min_seg_cm=ibdne_mincm)

#  remove ibd with tmrca < 1.5 (required according to IBDNe paper)
ibd.filter_ibd_by_time(min_tmrca=1.5)

# make diploid samples
# convert to heterzygous diploids
# Note: remove_hbd might not remove a lot segments as
# hdb only involves n/2 pairs of a total of n(n-1)/2 pairs (ratio: 1/(n-1)).
# When n is large, the difference is small
ibd.convert_to_heterozygous_diploids(remove_hbd=True)
ibd.flatten_diploid_ibd()

# calculate coverage and remove peaks
ibd.calc_ibd_cov()
ibd.find_peaks()

# ---- save ibd for ibdne call

# link ibdne.jar file
if not Path("ibdne.jar").exists():
    assert Path(args.ibdne_jar).exists()
    this = Path("ibdne.jar")
    target = Path(args.ibdne_jar).absolute()
    this.symlink_to(target)
    print(f"link {this} -> {target}")

# use nerunner (dry_run) to preare inputs and bash scripts
# for ibd before removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd,
    input_folder=".",
    output_folder=".",
    mincm=ibdne_mincm,
    minregion=ibdne_minregion,
)
nerunner.run(nthreads=6, mem_gb=20, dry_run=True)

# save of OBJ copy
ibd.pickle_dump(f"{genome_set_id}_orig.ibdne.ibdobj.gz")

# remove peaks

ibd2 = ibd.duplicate(f"{genome_set_id}_rmpeaks")
ibd2.remove_peaks()


# ---- save ibd for ibdne call
# for ibd after removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd2,
    input_folder=".",
    output_folder=".",
    mincm=ibdne_mincm,
    minregion=ibdne_minregion,
)
nerunner.run(nthreads=6, mem_gb=20, dry_run=True)

# save of OBJ copy
ibd.pickle_dump(f"{genome_set_id}_rmpeaks.ibdne.ibdobj.gz")


print(
    f"""
      output files:
        ibdne.jar
        {genome_set_id}_orig.sh
        {genome_set_id}_orig.map
        {genome_set_id}_orig.ibd.gz
        {genome_set_id}_orig.cov.pq
        {genome_set_id}_rmpeaks.sh
        {genome_set_id}_rmpeaks.map
        {genome_set_id}_rmpeaks.ibd.gz
        {genome_set_id}_ibddist.ibdobj.gz
        {genome_set_id}_orig_all.ibdcov.ibdobj.gz
        {genome_set_id}_orig_unrel.ibdcov.ibdobj.gz
        {genome_set_id}_orig.ibdne.ibdobj.gz
        {genome_set_id}_rmpeaks.ibdne.ibdobj.gz
    """
)
