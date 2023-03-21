#! /usr/bin/env python3
import argparse
import sys
import tsinfer
import tsdate
import msprime
import pandas as pd
import numpy as np
import allel
import io
from typing import Tuple

parser = argparse.ArgumentParser(
    """
    When one single-chromosome vcf file is provided, it infers treeseq per
    chromosome; when multiple single-chromome vcf files or single multiple-chromosme
    vcf file are provided, it infer treeseq genome-widely

    Output: {d.chrno}.trees (tree will have the same prefix as the vcf files)
    """
)
parser.add_argument("--vcf_files", nargs="+", type=str, required=True)
parser.add_argument("--r", type=float, required=True)
parser.add_argument("--u", type=float, required=True)
parser.add_argument(
    "--ne", type=str, required=True, help="interger or file name to true_ne table"
)

args = parser.parse_args()
print(args)
bp_per_cm = int(0.01 / args.r)
vcf_files = args.vcf_files
input_ne = args.ne

if input_ne.isnumeric():
    input_ne = int(input_ne)
else:
    # TODO: check true_ne table format for simulation code
    # Format: col1: generation, col2: Ne
    df_true_ne = pd.read_csv(input_ne, sep="\t")
    assert len(df_true_ne.columns) == 2
    df_true_ne.columns = ["Generation", "Ne"]
    assert pd.api.types.is_integer_dtype(df_true_ne.Generation)
    assert pd.api.types.is_integer_dtype(df_true_ne.Ne)
    # using harmonic mean for recent 100 generations
    input_ne = 1 / (1 / df_true_ne[df_true_ne.Generation <= 100].Ne).mean()
    input_ne = int(input_ne)  # tsdate complain if input_ne is numpy.int64


def read_vcf_files_into_df(vcf_files) -> Tuple[pd.DataFrame, pd.Series]:
    # ----------------- extract data from vcf files -----------------------
    gt, chrom, pos, ref, alt, contigs = [], [], [], [], [], []
    Samples = None
    for vcf in vcf_files:
        calldata = allel.read_vcf(vcf, alt_number=1)
        # assume vcf samples are haploid or pseudo homozygous diploid
        gt.append(calldata["calldata/GT"][:, :, 0])
        chrom.append(calldata["variants/CHROM"])
        pos.append(calldata["variants/POS"])
        ref.append(calldata["variants/REF"])
        alt.append(calldata["variants/ALT"])
        contigs.extend(
            [
                line
                for line in allel.read_vcf_headers(vcf).headers
                if line.startswith("##contig")
            ]
        )
        # If multiple vcf files are present, check sample names
        # are in the same order across these vcf files.
        if Samples is None:
            Samples = calldata["samples"]
        else:
            assert np.all(Samples == calldata["samples"])
    # ---------------- variant info ------------------------------------------
    df_variants = pd.DataFrame(
        {
            "CHROM": np.hstack(chrom),
            "POS": np.hstack(pos),
            "REF": np.hstack(ref),
            "ALT": np.hstack(alt),
        }
    )
    # ------------------- contig info ------------------------------------
    df_chrlen = (
        pd.read_csv(io.StringIO("".join(contigs)), header=None, sep="\t")
        .iloc[:, 0]
        .str.extract("##contig=<ID=(\w+),length=(\d+)>")
        # make sure length is interger, see sim_mut.py
    )
    df_chrlen.columns = ["CHROM", "Length"]
    # --------------- clean up and check consistency ---------------------
    df_chrlen = df_chrlen.sort_values(["CHROM"]).drop_duplicates(["CHROM"])
    unique_chrs = df_variants["CHROM"].unique().tolist()
    df_chrlen = df_chrlen[df_chrlen.CHROM.isin(unique_chrs)]
    print(df_chrlen, unique_chrs)
    assert df_chrlen.shape[0] == len(unique_chrs)
    df_chrlen["Length"] = df_chrlen.Length.astype(int)
    seqlens = df_chrlen.set_index("CHROM")["Length"]
    # ------------- add genome_wide position information---------------------
    chr_gw_starts = seqlens.cumsum() - seqlens
    chr_gw_starts = chr_gw_starts.rename("CHR_GW_START").reset_index()
    df_variants = df_variants.merge(chr_gw_starts, how="left", on="CHROM")
    df_variants["GW_POS"] = df_variants.POS + df_variants.CHR_GW_START
    # -------------- combine with GT data ----------------------------------
    df_variants = pd.DataFrame(
        np.vstack(gt), index=pd.MultiIndex.from_frame(df_variants)
    ).sort_index(level=["CHROM", "POS"])
    df_variants.columns = list(Samples)
    return df_variants, seqlens


df_variants, seqlens = read_vcf_files_into_df(vcf_files)

print(seqlens)

# ------- add samples to tsinfer --------------------------------
gw_seqlen = seqlens.sum()
with tsinfer.SampleData(
    sequence_length=gw_seqlen,
    path="tmp.genome_wide.samples",
    num_flush_threads=4,
) as ts_samples:
    Samples = list(df_variants.columns)
    # First add individual
    for sample_id in Samples:
        ts_samples.add_individual(
            ploidy=1, metadata={"name": sample_id}, population=None
        )
    # Then add site
    # here gt only use haplotype 0, due to pseudo-homogyzougs diploid
    Variants = df_variants.index.to_frame().reset_index(drop=True)[
        ["GW_POS", "REF", "ALT"]
    ]
    GT = df_variants.to_numpy()
    for (pos_site, ref_site, alt_site), gt_site in zip(
        Variants.itertuples(index=False), GT
    ):
        if pos_site > gw_seqlen:
            print(pos_site, ref_site, alt_site, gt_site[:10], file=sys.stderr)
            continue
        alleles = [ref_site, alt_site]
        ts_samples.add_site(
            position=pos_site,
            genotypes=gt_site,
            alleles=alleles,
        )

# ----- make rate map (with chr breaks) ------------------
r_chrom = 0.01 / bp_per_cm
r_break = np.log(2)
gw_starts = np.cumsum(seqlens) - seqlens
gw_ends = np.cumsum(seqlens)
gw_starts[1:] += 1
map_pos = np.sort(np.hstack([gw_starts, gw_ends])).tolist()
map_rate = [r_chrom, r_break] * len(seqlens)
map_rate = map_rate[:-1]
rate_map = msprime.RateMap(position=map_pos, rate=map_rate)


# --------- run tsinfer ---------------------------------
exclude_pos = None  # df_anc[df_anc.Reliable == 0].Pos.to_list()
sample_data = tsinfer.load("tmp.genome_wide.samples")
ts = tsinfer.infer(
    sample_data,
    recombination_rate=rate_map,
    mismatch_ratio=1 / 67,
    exclude_positions=exclude_pos,
)
ts.dump(f"tmp_genome_wide.trees")

# -------- Get the Ne parameter --------------------------


# --------- run tsdate ---------------------------------
ts_dated = tsdate.date(
    ts.simplify(keep_unary=False), Ne=input_ne, mutation_rate=args.u
)  # <============
ts_dated.dump("tmp.genome_wide_dated.trees")

# my tskibd.cpp do special thing when there is mutation information
# remove them for the purpose of IBD calling
num_sites = ts_dated.tables.sites.num_rows
ts_dated_sites_removed = ts_dated.delete_sites(np.arange(num_sites))
ts_dated_sites_removed.dump("tmp.genome_wide_dated_nosites.trees")

# --------- split GW tree to chromosome trees ------------------
df_chroms = pd.DataFrame(
    {
        "GW_START": seqlens.cumsum() - seqlens,
        "GW_END": seqlens.cumsum(),
    }
)
for chrom, chrom_start, chrom_end in df_chroms.itertuples():
    print(chrom_start, chrom_end, chrom_end - chrom_start)
    chrom_ts = ts_dated_sites_removed.keep_intervals(
        [[chrom_start, chrom_end]], simplify=True
    ).trim()
    chrom_ts.dump(f"{chrom}.trees")
