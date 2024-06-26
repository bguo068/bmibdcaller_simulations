#! /usr/bin/env python3

import ibdutils.utils.ibdutils as ibdutils
import ibdutils.runner.ibdne as ibdne
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from pathlib import Path
from subprocess import run

ibd_jar_default = str(Path(__file__).parent / "ibdne.jar")

# parse arguments
pa = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
pa.add_argument("--ibd_files", type=str, nargs=14, required=True)
pa.add_argument("--ibd_files_true", type=str, nargs=14, default=[])
pa.add_argument("--vcf_files", type=str, nargs=14, required=True)
pa.add_argument("--genome_set_id", type=str, required=True)
pa.add_argument("--ibdne_jar", type=str, default=ibd_jar_default)
# TODO: add the option to nextflow process and set a the default in params of nextflow
pa.add_argument(
    "--peak_validate_meth", type=str, choices=["xirs", "ihs"], default="xirs"
)
pa.add_argument("--ibdne_mincm", type=float, default=2)
pa.add_argument("--ibdne_minregion", type=float, default=10)
pa.add_argument("--ibdne_flatmeth", type=str, default="none")
# This can be used to specify which chromosome to use, especially when selection
# is not established on some chromosomes and we only use those that are
# established to better observe selection effects
pa.add_argument("--ibdne_chrlist", type=int, nargs="+", default=None)
pa.add_argument(
    "--ibdne_no_diploid_conversion",
    type=str,
    choices=["true", "false"],
    default="false",
)
args = pa.parse_args()

print(args)

assert args.ibdne_jar is not None
assert args.ibdne_flatmeth in ["none", "merge", "keep_hap_1_only"]


idx = args.genome_set_id
ibdne_mincm = args.ibdne_mincm
ibdne_minregion = args.ibdne_minregion

# output prefix string
label_str = "_".join(
    [
        str(x)
        for x in [
            args.genome_set_id,
            args.ibdne_mincm,
            args.ibdne_minregion,
            args.ibdne_flatmeth,
        ]
    ]
)

# read ibd
genome_14_100 = ibdutils.Genome.get_genome("simu_14chr_100cm")
ibd = ibdutils.IBD(genome=genome_14_100, label=f"{label_str}_orig")
ibd.read_ibd(ibd_fn_lst=args.ibd_files)
ibd.calc_ibd_cov()
ibd.find_peaks()


# output files:
of_ibddist_obj = f"{label_str}.ibddist.ibdobj.gz"


# store combined IBD for IBD distribution analysis
ibd.pickle_dump(of_ibddist_obj)


# make a matrix before possile loading of filterred IBD
mat = ibd.make_ibd_matrix()

# if the true IBD passed in, we can run tskibd-ibd filter.
# We then read the filter IBD (for IBDNe) and replace the IBD dataframe of the
# IBD object
ibd_files_filtered = []
if len(args.ibd_files_true) != 0:
    assert len(args.ibd_files_true) == 14
    # call tskibd-filter per chromosome
    for i, chrno in enumerate(range(1, 15)):
        inf = args.ibd_files[i]
        tru = args.ibd_files_true[i]
        out = f"filt_{chrno}.ibd"
        cmd = f"tskibd-filter {inf} {tru} -o {out} "
        run(cmd, shell=True, check=True)
        ibd_files_filtered.append(out)
    # load filterred IBD
    ibd_filt = ibdutils.IBD(genome=genome_14_100, label=f"{label_str}_orig")
    ibd_filt.read_ibd(ibd_fn_lst=ibd_files_filtered)
    # replace `_df` member of ibd variable but not the `coverage` and `peaks`
    ibd._df = ibd_filt._df

# remove highly relatedness samples
unrelated_samples = ibd.get_unrelated_samples(ibdmat=mat)
ibd.subset_ibd_by_samples(subset_samples=unrelated_samples)

# ######################################
# prepare input for IBDNe

#  remove ibd with tmrca < 1.5 (required according to IBDNe paper)
# this is only useful for tskibd ibd. IBD from other callers have tmrca set to 100
ibd.filter_ibd_by_time(min_tmrca=1.5)

ibd.filter_ibd_by_length(min_seg_cm=ibdne_mincm)

# calculate XiR,s
if args.peak_validate_meth == "xirs":
    ibd.calc_xirs(vcf_fn_lst=args.vcf_files, min_maf=0.01)
elif args.peak_validate_meth == "ihs":
    ibd.calc_ihs(vcf_fn_lst=args.vcf_files, min_maf=0.01)
else:
    raise NotImplementedError(f"{args.peak_validate_meth} is not valid method")


if args.ibdne_no_diploid_conversion == "true":
    print("skip diploid conversion")
elif args.ibdne_no_diploid_conversion == "false":
    # convert to heterzygous diploids
    # Note: remove_hbd might not remove a lot segments as
    # hdb only involves n/2 pairs of a total of n(n-1)/2 pairs (ratio: 1/(n-1)).
    # When n is large, the difference is small
    ibd.convert_to_heterozygous_diploids(remove_hbd=True)

    if args.ibdne_flatmeth != "none":
        ibd.flatten_diploid_ibd(method=args.ibdne_flatmeth)
else:
    raise NotImplementedError(
        f"{args.ibdne_no_diploid_conversion} is not valid value is not a value value for option --ibdne_no_diploid_conversion"
    )

# calculate coverage and remove peaks
ibd.calc_ibd_cov()
ibd.find_peaks()


if args.peak_validate_meth == "xirs":
    xirs_df = ibd._xirs_df
    ibd.filter_peaks_by_xirs(xirs_df)
elif args.peak_validate_meth == "ihs":
    ibd.filter_peaks_by_ihs(min_ihs_hits=1)
else:
    raise NotImplementedError(f"{args.peak_validate_meth} is not valid method")

ibd2 = ibd.duplicate(f"{label_str}_rmpeaks")
ibd2.remove_peaks()
ibd2._df = ibd2.cut_and_split_ibd()


# filter by chromosome for Ne
if args.ibdne_chrlist is not None:
    ibd._df = ibd._df[lambda df: df.Chromosome.isin(args.ibdne_chrlist)]
    ibd2._df = ibd2._df[lambda df: df.Chromosome.isin(args.ibdne_chrlist)]

# write ne file
of_orig_ibdne_obj = f"{label_str}_orig.ibdne.ibdobj.gz"
ibd.pickle_dump(of_orig_ibdne_obj)

of_rmpeaks_ibdne_obj = f"{label_str}_rmpeaks.ibdne.ibdobj.gz"
ibd2.pickle_dump(of_rmpeaks_ibdne_obj)

# link ibdne.jar file
if not Path("ibdne.jar").exists():
    assert Path(args.ibdne_jar).exists()
    this = Path("ibdne.jar")
    target = Path(args.ibdne_jar).absolute()
    this.symlink_to(target)
    print(f"link {this} -> {target}")

# use nerunner (dry_run) to preare inputs and bash scripts
# --- for ibd before removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd,
    input_folder=".",
    output_folder=".",
    mincm=ibdne_mincm,
    minregion=ibdne_minregion,
)
nerunner.run(nthreads=10, mem_gb=20, dry_run=True)

# --- for ibd after removing ibd peaks ------------------------
nerunner = ibdne.IbdNeRunner(
    ibd=ibd2,
    input_folder=".",
    output_folder=".",
    mincm=ibdne_mincm,
    minregion=ibdne_minregion,
)
nerunner.run(nthreads=10, mem_gb=20, dry_run=True)

print(
    f"""
      output files:
        {of_ibddist_obj}
        ibdne.jar
        {idx}_orig.sh
        {idx}_orig.map
        {idx}_orig.ibd.gz
        {idx}_rmpeaks.sh
        {idx}_rmpeaks.map
        {idx}_rmpeaks.ibd.gz
        {of_orig_ibdne_obj}
        {of_rmpeaks_ibdne_obj}
      """
)
