# -*- coding: utf-8 -*-
"""
Created: 08-04-2023
Author Anna Schedl
Description: This script takes an output file from "virpool_v763c212" and converts it into a
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
        python script.py -virpool_v763c212 -virpool_v763c212_summary.csv -min_threshold
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

virpool_df = pd.DataFrame({"timepoint": pd.Series(dtype=int)})
for f in args.input:
        with open(f, 'r') as file:
            df = pd.json_normalize(yaml.safe_load(file))
            # timepoint is found in file name after "-"
            s = f[f.rfind("/")+1:f.rfind(".")]
            df["timepoint"] = s[s.find("-")+1:]
            df["timepoint"] = df["timepoint"].astype(int)
            df.rename(columns={"other": "others"}, inplace=True)
        virpool_df = pd.concat([virpool_df, df], ignore_index=True).sort_values(by=["timepoint"])

if virpool_df.empty:
      with open(args.output, 'w') as f:
            pass
else:
      # filter dataset
      virpool_df = func.filter_dataframe(virpool_df.copy(), args.minThreshold, tool_name=args.name.lower())

      # write to summary file
      virpool_df.to_csv(args.output, index=False)