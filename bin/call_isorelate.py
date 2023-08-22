#! /usr/bin/env python3
import argparse
from pathlib import Path
from subprocess import run


def get_args(cmd_args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", type=str, required=True)
    parser.add_argument("--r", type=float, required=True)
    parser.add_argument("--seqlen", type=int, required=True)
    parser.add_argument("--chrno", type=int, required=True)
    parser.add_argument("--min_snp", type=int, default=20)
    parser.add_argument("--min_len_bp", type=int, default=50000)
    parser.add_argument("--minmac", type=int, required=True)
    parser.add_argument("--imiss", type=float, required=True)
    parser.add_argument("--vmiss", type=float, required=True)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--genome_set_id", type=int, required=True)
    if cmd_args is not None:
        return parser.parse_args(cmd_args)
    else:
        return parser.parse_args()


args = get_args()
vcf = args.vcf
cpus = args.cpus
min_snp = args.min_snp
min_len_bp = args.min_len_bp
bp_per_cm = 0.01 / args.r
seqlen_in_cm = args.seqlen / bp_per_cm
# use `minmac` instead of `maf` to make it consistent across different IBD callers
minmac = args.minmac
imiss = args.imiss
vmiss = args.vmiss
chrno = args.chrno

# convert vcf to ped map via plink and modify the cM in the map file
cmd = f"""
plink --vcf {vcf} --recode ped 12  -out ped --double-id
mv ped.map pre_fix_ped.map
awk -v bp_per_cm={bp_per_cm} '{{$3=$4/bp_per_cm; print}}' pre_fix_ped.map > ped.map
rm pre_fix_ped.map
"""
assert run(cmd, shell=True).returncode == 0

ofn = f"{args.genome_set_id}_{chrno}_isorelate.ibd"

# NOTE: treat all sample as haploid, set 5th col of ped file to 1
r_code = f"""
ped <- read.delim("ped.ped", sep=' ', header=F)
ped[, 5] <- 1
nsample = nrow(ped)
maf = {minmac} / 2.0 / nsample
map <- read.delim("ped.map", sep=' ', header=F)
ped_map <- list(ped=ped, map=map)
library(isoRelate)
# reformat and filter genotypes
my_genotypes <- getGenotypes(ped.map = ped_map,
                             reference.ped.map = NULL,
                             maf = maf,
                             isolate.max.missing = {imiss},
                             snp.max.missing = {vmiss},
                             chromosomes = NULL,
                             input.map.distance = "cM",
                             reference.map.distance = "cM")
# estimate parameters
my_parameters <- getIBDparameters(ped.genotypes = my_genotypes,
                                   number.cores = {cpus})
my_ibd <- getIBDsegments(ped.genotypes = my_genotypes,
                         parameters = my_parameters,
                         number.cores = {cpus},
                         minimum.snps = {min_snp},
                         minimum.length.bp = {min_len_bp},
                         error = 0.001)
# reformat
my_ibd <- my_ibd[c('iid1', 'iid2', 'start_position_bp', 'end_position_bp')]
colnames(my_ibd) <- c("Id1", "Id2", "Start", "End")
my_ibd$Id1 = gsub('tsk_', '', my_ibd$Id1)
my_ibd$Id2 = gsub('tsk_', '', my_ibd$Id2)
my_ibd$Ancestor = 9999
my_ibd$Tmrca = 100
my_ibd$HasMutation = 0
write.table(my_ibd, "{ofn}", sep='\t', quote=F, row.names=F)
"""
Path("run_isorelate.R").write_text(r_code)
# ------- run R code to generate the ibd file
cmd = ["Rscript", "run_isorelate.R"]
assert run(cmd).returncode == 0

# # ------------------- make genetic map ----------------------------
with open(f"{chrno}.map", "w") as f:
    fields = [f"{chrno}", ".", "0", "1"]
    f.write("\t".join(fields) + "\n")
    fields = [f"{chrno}", ".", f"{seqlen_in_cm}", f"{bp_per_cm * seqlen_in_cm}"]
    f.write("\t".join(fields) + "\n")

# -------------------- write log (just for consistency with other ibd callers) -
Path(f"{chrno}.log").write_text(" ".join(cmd))
