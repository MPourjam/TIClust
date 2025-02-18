import pathlib as pl

config_path = pl.Path(__file__).resolve().parent.parent.joinpath("bin/vsearch")

tic_configs = {
    "VSEARCH_BIN_PATH": str(config_path)
}
