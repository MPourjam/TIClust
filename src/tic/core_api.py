import argparse
from pathlib import Path
from os import cpu_count


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta-file",
        help="Fasta file with taxonomies starting with 'tax='"
    )
    parser.add_argument(
        "-z",
        "--zotu-table",
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
        "-ft"
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
    args = parser.parse_args()
    input_fasta_file = Path(args.fasta_file).resolve()
    zotu_table = Path(args.zotu_table).resolve()
    args.fasta_file = str(input_fasta_file)
    args.zotu_table = str(zotu_table)

    return args
