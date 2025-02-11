import sys
import argparse
import pathlib as pl
from .simple_tic import TICAnalysis
from .core_api import parse_arguments

def main():
    args = parse_arguments()
    simple_tic = TICAnalysis(args.fasta_file)
    simple_tic.run()

if __name__ == "__main__":
    main()
