#! /usr/bin/env python3
import argparse
from pathlib import Path
import sys
import tskit
import msprime
import pandas as pd
import gzip

# usage:
#    sim_mut.py \
#         --in_chrno ${d.chrno} \
#         --in_trees ${in_trees_fn} \
#         --mu ${d.u} \
#         --out_vcf ${d.chrno}.vcf.gz

parser = argparse.ArgumentParser()
parser.add_argument("--in_chrno", type=int, required=True)
parser.add_argument("--in_trees", type=str, required=True)
parser.add_argument("--mu", type=float, required=True)
parser.add_argument("--out_vcf", type=str, required=True, help="xx.vcf.gz")

args = parser.parse_args()
chrno = args.in_chrno
trees_fn = args.in_trees
mu = args.mu
out_vcf = args.out_vcf

# load trees
ts = tskit.load(trees_fn)

# # for testing:
# ts = msprime.sim_ancestry(
#     100,
#     population_size=1000,
#     recombination_rate=0.01 / 15000,
#     sequence_length=100 * 15000,
# )
# ts.dump("tmp.trees")

# remove exising slim mutation (it has different ref/alt coding)
site_list = list(range(ts.num_sites))
ts = ts.delete_sites(site_list)

# add neutral mutation
ts_mutated = msprime.sim_mutations(ts, rate=mu, random_seed=chrno)

# write pseudo homozygous vcf files


def write_peudo_homozygous_vcf(ts_mutated, out_vcf):
    # extract information from ts
    gt_list = []
    pos_list = []
    ref_list = []
    alt_list = []
    for v in ts_mutated.variants():
        if len(v.alleles) != 2:
            continue
        gt_list.append(v.genotypes)
        ref_list.append(v.alleles[0])
        alt_list.append(v.alleles[1])
        pos_list.append(int(v.position))

    # prep header
    header = f"""##fileformat=VCFv4.2
##source=tskit 0.4.0
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID={chrno},length={int(ts_mutated.sequence_length)}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""
    # colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    # variant info
    df1 = pd.DataFrame(
        {
            "#CHROM": chrno,
            "POS": pos_list,
            "ID": ".",
            "REF": ref_list,
            "ALT": alt_list,
            "QUAL": ".",
            "FILTER": "PASS",
            "INFO": ".",
            "FORMAT": "GT",
        }
    )

    # gt info
    df2 = pd.DataFrame(gt_list)
    df2.columns = [f"tsk_{n}" for n in df2.columns]
    df2 = df2.astype(str)
    df2 = df2 + "|" + df2

    # combine variant info with gt infor
    df = pd.concat([df1.reset_index(drop=True), df2.reset_index(drop=True)], axis=1)

    # write
    with gzip.open(out_vcf, "wt") as f:
        f.write(header)
        df.to_csv(f, sep="\t", header=True, index=False)


write_peudo_homozygous_vcf(ts_mutated, out_vcf)
