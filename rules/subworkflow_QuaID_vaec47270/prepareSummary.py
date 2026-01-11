# -*- coding: utf-8 -*-
"""
Created: 08-04-2024
Author Anna Schedl
Description: This script takes an output file from "QuaID_vaec47270" and converts it into a
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
        python script.py -QuaID_vaec47270 -QuaID_vaec47270_summary.csv -min_threshold
"""
import sys, argparse, yaml
import pandas as pd

sys.path.append("bin/custom_scripts")
import PostPred_functions as func  # pyright: ignore[reportMissingImports]

# create Parser object
parser = argparse.ArgumentParser()

# Add arguments
parser.add_argument("-i", "--input", nargs="+", help="List of input file names")
parser.add_argument("--metaFile", type=str, help="meta data file")
parser.add_argument("--whoFile", type=str, help="WHO VOC definition file in yaml format.")
parser.add_argument("--name", type=str, help="Tool Name")
parser.add_argument("--minThreshold", type=float, help="Minimal threshold for relative abundance")
parser.add_argument("-o", "--output", type=str, help="output file name")

args = parser.parse_args()

# read out data
meta_df = pd.read_csv(args.metaFile, sep='\t')
meta_df["sample_date"] = pd.to_datetime(meta_df["sample_date"], format='%Y-%m-%d')

data = []
for f in args.input:
        data.append(pd.read_csv(f, sep=','))

quaid_df = pd.concat(data)

if quaid_df.empty: # create empty file in case result files were empty
        with open(args.output, 'w') as f: pass

else:
        quaid_df["sample_date"] = pd.to_datetime(quaid_df["Date"], format='%Y-%m-%d')
        
        quaid_df = quaid_df[quaid_df["Total AF (quasi-unique)"] != 0.0]
        quaid_df = quaid_df.merge(meta_df, on="sample_date")

        # translate WHO names to primary pangolin lineage associated with it
        quaid_df["pango"] = ""

        for index, row in quaid_df.iterrows():
                pango = func.who_pango_translate(args.whoFile, row["WHO name"])

                if len(pango) != 0:
                    quaid_df.at[index, "pango"] = pango
                else:
                    quaid_df.at[index, "pango"] = row["WHO name"]
        
        # change Total AF to percentage range [0,1]
        quaid_df['Total AF'] = quaid_df['Total AF (quasi-unique)'] / quaid_df['Total QU count']

        quaid_df = quaid_df.pivot(index='timepoint', columns='pango', values='Total AF')
        quaid_df.reset_index(inplace=True)
        quaid_df.columns.name = None

        # filter dataset
        quaid_df = func.filter_dataframe(quaid_df.copy(), args.minThreshold, tool_name=args.name.lower())

        # write to summary file
        quaid_df.to_csv(args.output, index=False)