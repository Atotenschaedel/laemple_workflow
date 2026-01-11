# -*- coding: utf-8 -*-
"""
Created: 29-06-2023
Author Anna Schedl
Description: This script takes an output file from "vaquero_v93e7db3" and converts it into a
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
        python script.py -vaquero_v93e7db3 -vaquero_v93e7db3_summary.csv -min_threshold
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

vaquero_df = pd.DataFrame({'timepoint': pd.Series(dtype=int)})
meta_df = pd.read_csv(args.metafile, sep='\t')
meta_df["sample_date"] = pd.to_datetime(meta_df["sample_date"], format='%Y-%m-%d')

with open(args.input, 'r') as f:
        first_char = f.read(1)
        if not first_char:
                # write emtpy output file if no results were created
                with open(args.output, 'w') as f:
                        pass
        else:
                vaquero_df = pd.read_csv(args.input, sep='\t')
                vaquero_df["sample_date"] = pd.to_datetime(vaquero_df["sample_date"], format='%Y-%m-%d')
                vaquero_df = vaquero_df[vaquero_df["value"] != 0.0]
                vaquero_df = vaquero_df.merge(meta_df, on='sample_date')
                vaquero_df.drop(columns=["sample_id", "sample", "sample_date", "LocationID", "LocationName"], inplace=True)
                vaquero_df = vaquero_df.pivot(index='timepoint',columns='variant', values='value')
                vaquero_df.reset_index(inplace=True)
                vaquero_df.columns.name = None

                # filter dataset
                vaquero_df = func.filter_dataframe(vaquero_df.copy(), args.minThreshold, tool_name=args.name.lower())

                # write to summary file
                vaquero_df.to_csv(args.output, index=False)