#!/usr/bin/env python

import argparse
import os
import sys
import re
from pathlib import Path

import pandas as pd


def parse() -> tuple:
    parser = argparse.ArgumentParser(
        description="Demultiplex .fastq files given the sample and experiment barcodes."
    )
    parser.add_argument("-b", type=str, help="Barcode table (csv) path")
    parser.add_argument("-f", type=str, help=".fastq file")

    _args = parser.parse_args()
    return _args.b, _args.f

def validate_inputs(bc: str, fs: str) -> tuple:
    if not (bc_path := Path(bc)).exists():
        sys.stderr.write("InputError: Barcode .csv does not exist\n")
        quit()

    if not (file_path := Path(fs)).exists():
        sys.stderr.write("InputError: .fastq file does not exist\n")
        quit()

    bc_df = pd.read_csv(bc_path, sep=";").dropna(axis=0, inplace=False)
    bc_df.Library.replace(" ", "_", inplace=True, regex=True)

    return bc_df, file_path

def demultiplex(bc_csv: Path, f1: Path) -> None:
    bcs: list = [
        f"{(row[1]['Library']).strip()}={(row[1]['Sequence']).strip()}$"
        for row in bc_csv.iterrows()
    ]
    bcs = " -a ".join(bcs)  # 3' anchored i7 barcodes
    
    # input file name example = "forward_R1.fastq.gz"
    # pref = f1.name.split(".")[0].split("_")[0] # forward/reverse 

    cmd = f"""cutadapt -e 0.2 --no-indels -j 0 -a {bcs} -o {{name}}.fastq.gz {f1}"""
    os.system(cmd)


if __name__ == "__main__":
    args = parse()
    bc_df, file = validate_inputs(*args)
    demultiplex(bc_df, file)
