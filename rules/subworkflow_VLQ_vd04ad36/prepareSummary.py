# -*- coding: utf-8 -*-
"""
Created: 29-05-2023
Author Anna Schedl
Description: This script takes an output file from "VLQ_vd04ad36" and converts it into a
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
        python script.py -VLQ_vd04ad36 -VLQ_vd04ad36_summary.csv -min_threshold
"""
import sys, argparse, yaml
import pandas as pd

sys.path.append("bin/custom_scripts")
import PostPred_functions as func  # pyright: ignore[reportMissingImports]

# create Parser object
parser = argparse.ArgumentParser()

# Add arguments
parser.add_argument("-i", "--input", nargs="+", help="List of input file names")
parser.add_argument("--name", type=str, help="Tool Name")
parser.add_argument("--minThreshold", type=float, help="Minimal threshold for relative abundance")
parser.add_argument("-o", "--output", type=str, help="output file name")

args = parser.parse_args()

# read out result files
VLQ_df = pd.DataFrame({'timepoint': pd.Series(dtype=int)})
for f in args.input:
        #colum names: [variant, tpm, freq(%), adj_freq(%)]
        df = pd.read_csv(f, sep='\t', skiprows=3, header=None, names=["variant", "tpm", "freq(%)", "adj_freq(%)"])
        df = df[df["adj_freq(%)"] != 0.0]
        df.reset_index(inplace=True, drop=True)
        df.drop(columns=["tpm", "freq(%)"], inplace=True)
                    
        # adjust to abundance between 1-0
        df["adj_freq(%)"] = df["adj_freq(%)"] / 100

        df = df.T
        df.columns = df.iloc[0]
        df = df[1:]
        df.reset_index(inplace=True, drop=True)

        # timepoint is found in file name after "-"
        s = f[f.rfind("/")+1:f.rfind(".")]
        df["timepoint"] = s[s.find("-")+1:s.rfind("_")]

        VLQ_df = pd.concat([VLQ_df, df], ignore_index=True).sort_values(by=["timepoint"])

        
if VLQ_df.empty: # create empty file in case result files were empty
        with open(args.output, 'w') as f:
              pass
else:
        # filter dataset
        VLQ_df = func.filter_dataframe(VLQ_df.copy(), args.minThreshold , tool_name=args.name.lower())

        # write to summary file
        VLQ_df.to_csv(args.output, index=False)