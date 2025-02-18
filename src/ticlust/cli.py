from .ticlust import TICAnalysis
from .core_api import parse_arguments

def main():
    args = parse_arguments()
    simple_tic = TICAnalysis(
        args.fasta_file,
        args.zotu_table
        )
    threshold = {
        "species": args.species_thr,
        "genus": args.genera_thr,
        "family": args.family_thr
    }
    simple_tic.run(
        threads = args.threads,
        cluster_thresholds_d = threshold
    )

if __name__ == "__main__":
    main()
