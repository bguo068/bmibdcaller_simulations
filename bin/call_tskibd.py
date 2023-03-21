#! /usr/bin/env python3
from subprocess import run
import pandas as pd
import argparse


def get_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--tree", type=str, required=True)
    parser.add_argument("--chrno", type=int, required=True)
    parser.add_argument("--r", type=float, default=0.01 / 15000)
    parser.add_argument("--seqlen", type=int, default=15000 * 100)
    parser.add_argument("--mincm", type=float, default=2.0)
    parser.add_argument("--genome_set_id", type=int, required=True)
    parser.add_argument(
        "--max_tmrca", type=int, nargs="*", help="more than 1 cause filtering"
    )
    args = parser.parse_args()

    args.bp_per_cm = int(0.01 / args.r)
    args.seqlen_in_cm = args.seqlen * 100 * args.r
    print(args)

    return args


# ------------------- make genetic map ----------------------------
def write_map_file(ofs_map, args):
    with open(ofs_map, "w") as f:
        fields = [f"{args.chrno}", ".", "0", "1"]
        f.write("\t".join(fields) + "\n")
        fields = [
            f"{args.chrno}",
            ".",
            f"{args.seqlen_in_cm}",
            f"{args.seqlen}",
        ]
        f.write("\t".join(fields) + "\n")


def run_tskibd(ofs_tskibd, args):
    sample_window = int(0.01 * args.bp_per_cm)
    assert (
        run(
            # run tskibd within subfolder to avoid folder containmination
            f"""
                mkdir tskibd_{args.chrno}
                cd tskibd_{args.chrno}
                tskibd {args.chrno} {args.bp_per_cm} {sample_window} {args.mincm} \
                    ../{args.tree}
                cd ..
                mv tskibd_{args.chrno}/{args.chrno}.ibd {ofs_tskibd}
                rm -rf tskibd_{args.chrno}
            """,
            shell=True,
            check=True,
        ).returncode
        == 0
    )


if __name__ == "__main__":
    args = get_args()

    # out file names
    ofs_map = f"{args.genome_set_id}_{args.chrno}.map"
    ofs_tskibd = f"{args.genome_set_id}_{args.chrno}_tskibd.ibd"

    write_map_file(ofs_map, args)

    run_tskibd(ofs_tskibd, args)

    extra_files = []
    if args.max_tmrca is not None:
        for m in args.max_tmrca:
            ofn = f"{args.genome_set_id}_{args.chrno}_tskibd.ibd_maxtmrca_{m}"
            df = pd.read_csv(ofs_tskibd, sep="\t")
            df[df.Tmrca <= m].to_csv(ofn, index=None, sep="\t")
            extra_files.append(ofn)

    print(
        f"""
        output files:
            {ofs_tskibd}
            {ofs_map} 

        extra files:
            {extra_files}
            
        """
    )
