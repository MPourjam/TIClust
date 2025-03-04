import argparse
from pathlib import Path
from os import cpu_count
from ticlust import __version__

VERSION = __version__


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta-file",
        help="Fasta file with taxonomies starting with 'tax='",
        required=True
    )
    parser.add_argument(
        "-z",
        "--zotu-table",
        default="",
        help="zOTUs table file path"
    )
    parser.add_argument(
        "-st",
        "--species-thr",
        default=0.987,
        help="Similarity threshold for species-level clustering"
    )
    parser.add_argument(
        "-gt",
        "--genera-thr",
        default=0.95,
        help="Similarity threshold for genus-level clustering"
    )
    parser.add_argument(
        "-ft",
        "--family-thr",
        default=0.9,
        help="Similarity threshold for family-level clustering"
    )
    parser.add_argument(
        "-t",
        "--threads",
        default=int(cpu_count() * 0.75),
        help="Number of threads to use"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}"
    )
    args = parser.parse_args()
    input_fasta_file = Path(args.fasta_file).resolve()
    zotu_table = str(Path(args.zotu_table).resolve()) if args.zotu_table else None
    args.fasta_file = str(input_fasta_file)
    args.zotu_table = str(zotu_table)
    args.threads = int(args.threads)
    args.species_thr = float(args.species_thr)
    args.genera_thr = float(args.genera_thr)
    args.family_thr = float(args.family_thr)

    return args
