#!/usr/bin/env python

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


def parse() -> tuple:
    parser = argparse.ArgumentParser(
        description="Demultiplex .fastq files given the sample and experiment barcodes."
    )
    parser.add_argument("-b", type=str, help="Barcode table (csv) path")
    parser.add_argument("-f", type=str, help=".fastq files directory")

    _args = parser.parse_args()
    return _args.b, _args.f

def validate_inputs(bc: str, fs: str) -> tuple:
    if not (bc_path := Path(bc)).exists():
        sys.stderr.write("InputError: Barcode csv does not exist\n")
        quit()

    if not (fs_path := Path(fs)).exists():
        sys.stderr.write("InputError: .fastq files directory does not exist\n")
        quit()

    bc_df = pd.read_csv(bc_path, sep=";").dropna(axis=0, inplace=False)
    bc_df.Sample.replace(" ", "_", inplace=True, regex=True)

    fastq = next(fs_path.glob("*.fastq.gz"))

    return bc_df, fastq

def demultiplex(bc_csv: Path, fastq: Path) -> None:
    # using Reverse barcode for fw reads
    bcs: list = [
        f"{(row[1]['Sample']).strip()}={(row[1]['Reverse']).strip()}$"
        for row in bc_csv.iterrows()
    ]
    bcs = " -a ".join(bcs)  # anchored 3' adapters

    # input example = "trimmed_forward_libN.fastq.gz" 
    prefix = fastq.name.split(".")[0].removeprefix("trimmed_") # forward_libN

    cmd = f"""cutadapt -e 0.2 --no-indels -j 0 -a {bcs} -o {prefix}_{{name}}.fastq.gz {fastq}"""
    os.system(cmd)


if __name__ == "__main__":
    args = parse()
    bc_df, f_path = validate_inputs(*args)
    demultiplex(bc_df, f_path)
