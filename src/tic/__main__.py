from .simple_tic import TICAnalysis
from .core_api import parse_arguments


def run_simple_tic(input_fasta_file):
    tic_anlysis = TICAnalysis(input_fasta_file)
    tic_anlysis.run()

def main():
    args = parse_arguments()
    run_simple_tic(args.fasta_file)
