# -*- coding: utf-8 -*-
"""
Created: 29-05-2023
Author Anna Schedl
Description: This script takes an output file from "freyja_v1.4.3" and converts it into a
    standardized output format.

    The standardized format includes the following columns:
        - timepoint
        - lineage 1
        - lineage 2
        - ...
        - others
        - sample_name
        - tool_name

Example usage:
        python script.py -freyja_v1.4.3 -freyja_summary.csv -min_threshold
"""
import sys, argparse
import pandas as pd


sys.path.append("bin/custom_scripts")
import PostPred_functions as func  # pyright: ignore[reportMissingImports]


# create Parser object
parser = argparse.ArgumentParser()

# Add arguments
parser.add_argument("-i", "--input", nargs="+", help="List of input file names")
parser.add_argument("-n", "--name", help="Tool Name")
parser.add_argument("-m", "--minThreshold", help="Minimal threshold for relative abundance")
parser.add_argument("-o", "--output", help="output file name")

args = parser.parse_args()


frey_df = pd.DataFrame({"timepoint": pd.Series(dtype=int)})

for f in args.input:
    names = []
    values = []
    with open(f) as file:
        lines = file.readlines()
        for line in lines:
            line = line.replace('\n', "")
            if line.startswith("lineages"):
                line = line.replace("lineages\t", "")
                names.extend(line.split(" "))
            elif line.startswith("abundances"):
                line = line.replace("abundances\t", "")
                values.extend(line.split(" "))
    
    df = pd.DataFrame([values], columns=names)

    # Timepoint is number between last "-" and ".tsv"
    df["timepoint"] = f[f.rfind("-")+1:f.find(".tsv")]
    df["timepoint"] = df["timepoint"].astype(int)

    frey_df = pd.concat([frey_df, df], ignore_index=True).sort_values(by=["timepoint"])

if df.empty:
    # write empty output file if no results were created
    with open(args.output, 'w') as f:
        pass
else:
    # filter dataset
    frey_df = func.filter_dataframe(frey_df.copy(), float(args.minThreshold), tool_name=args.name)

    # write to summary file
    frey_df.to_csv(args.output, index=False)