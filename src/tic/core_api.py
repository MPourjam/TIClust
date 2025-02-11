import argparse
from pathlib import Path
from .simple_tic import TICAnalysis


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta-file",
        help="Fasta file with taxonomies starting with 'tax='"
    )
    args = parser.parse_args()
    input_fasta_file = Path(args.fasta_file).resolve()
    args.fasta_file = str(input_fasta_file)

    return args
