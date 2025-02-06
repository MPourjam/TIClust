import sys
import argparse
import pathlib as pl
from simple_tic import simple_tic as stic

def main():
    if len(sys.argv) < 2:
        print("Usage: simpletic <arg1>")
        sys.exit(1)
    
    # -f, --fasta-file
    # -o, --output-file
    # -t, --threshold
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta-file", help="Fasta file")
    # parser.add_argument("-o", "--output-file", help="Output file")
    # parser.add_argument("-t", "--threshold", help="Threshold")

    args = parser.parse_args()
    input_fasta_file = pl.Path(args.fasta_file).absolute().resolve()
    simple_tic = stic.TICAnalysis(input_fasta_file)
    simple_tic.run()

if __name__ == "__main__":
    main()