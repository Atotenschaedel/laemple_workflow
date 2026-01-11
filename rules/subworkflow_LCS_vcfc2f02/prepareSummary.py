# -*- coding: utf-8 -*-
"""
Created: 27-08-2024
Author Anna Schedl
Description: This script takes an output file from "LCS_vcfc2f02" and converts it into a
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
        python script.py -LCS_vcfc2f02 -LCS_vcfc2f02_summary.csv -min_threshold
"""
import sys, argparse, yaml
import pandas as pd

sys.path.append("bin/custom_scripts")
import PostPred_functions as func  # pyright: ignore[reportMissingImports]

# create Parser object
parser = argparse.ArgumentParser()

# Add arguments
parser.add_argument("-i", "--input", help="Input File")
parser.add_argument("--name", type=str, help="Tool Name")
parser.add_argument("--minThreshold", type=float, help="Minimal threshold for relative abundance")
parser.add_argument("-o", "--output", type=str, help="output file name")

args = parser.parse_args()



lcs_df = pd.read_csv(args.input, sep='\t')

if lcs_df.empty: # create empty file in case result files were empty
        with open(args.output, 'w') as f:
              pass

else:
        lcs_df[['exp', 'timepoint']] = lcs_df['sample'].str.split("-", expand = True)

        for index, row in lcs_df.iterrows():
                if '_' in row['variant_group']: lcs_df[index, 'variant_group'] = row['variant_group'][(row['variant_group'].find("_"))+1:]
                if '_' in row['variant_group']: lcs_df[index, 'variant_group'] = row['variant_group'][:row['variant_group'].find("_")]
                
        lcs_df = lcs_df.pivot(index='timepoint', columns='variant_group', values='proportion')
        lcs_df.reset_index(inplace=True)
        lcs_df.columns.name = None
        lcs_df['timepoint'] = lcs_df['timepoint'].astype(int)

        rename_dict = {}
        for col in lcs_df.columns:
                if "_" in col: 
                        var = col[col.find("_")+1:]
                        if "_" in var: rename_dict[col] = var[:var.find("_")]
                        else: rename_dict[col] = var

        lcs_df.rename(columns=rename_dict, inplace=True)

        # filter dataset
        lcs_df = func.filter_dataframe(lcs_df.copy(), args.minThreshold, tool_name=args.name.lower())

        # write to summary file
        lcs_df.to_csv(args.output, index=False)