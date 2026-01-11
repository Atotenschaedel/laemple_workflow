# -*- coding: utf-8 -*-
"""
Created: 29-05-2023
Author Anna Schedl
Description: This script takes an output file from "lollipop_v0.3.0" and converts it into a
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
        python script.py -lollipop_v0.3.0 -lollipop_v0.3.0_summary.csv -min_threshold
"""
import sys, argparse
import pandas as pd

sys.path.append("bin/custom_scripts")
import PostPred_functions as func  # pyright: ignore[reportMissingImports]


# create Parser object
parser = argparse.ArgumentParser()

# Add arguments
parser.add_argument("-i", "--input", type=str, help="List of input file names")
parser.add_argument("--metafile", type=str, help="meta data file")
parser.add_argument("--name", type=str, help="Tool Name")
parser.add_argument("--minThreshold", type=float, help="Minimal threshold for relative abundance")
parser.add_argument("-o", "--output", type=str, help="output file name")

args = parser.parse_args()


# read out input file 
df = pd.read_csv(args.input, sep='\t')
df["date"] = pd.to_datetime(df["date"], format='%Y-%m-%d')

# read out metadata file
meta_df = pd.read_csv(args.metafile, sep='\t')
meta_df["sample_date"] = pd.to_datetime(meta_df["sample_date"], format='%Y-%m-%d')

df = df.merge(meta_df, left_on ='date', right_on='sample_date')
df.drop(columns=["location", "date", "sample", "sample_date"], inplace=True)
df = df.pivot(index='timepoint',columns='variant', values='proportion')
df.rename(columns = {"undetermined": "others"}, inplace=True)
df.reset_index(inplace=True)
df.columns.name = None

if df.empty:
        # write empty output file if no results were created
        with open(args.output, 'w') as f:
                pass
else:
        # filter dataset
        df = func.filter_dataframe(df.copy(), args.minThreshold, tool_name=args.name.lower())

        # write to summary file
        df.to_csv(args.output, index=False)