# -*- coding: utf-8 -*-
"""
Created on 29-07-2023
Author: aschedl
Description: This script takes output from simulation and variant calling pipoeline
    and coverts it into a summary file.
    
    The output file includes the following columns:
        - timepoint
        - lineage 1 (relative abundance 0-1)
        - lineage 2 (relative abundance 0-1)
        - ...
        - others (relative abundance 0-1)
        - sample_name
        - tool_name 
        - sample
        - sample_date
        - coverage_avg (overall coverage average)
        - coverage_sd
        - uniformity_wg_per (uniformity overall - percentage of bases in all targeted regions (or whole‑genome))
        - MAPQ_avg (average quality score of mapped reads)
"""
import argparse
import pandas as pd
import numpy as np
import PostPred_functions as func

# create Parser object
parser = argparse.ArgumentParser()

# Add arguments
parser.add_argument("-s", "--samples", nargs="+", help="List of sample names")
parser.add_argument("-d", "--data_file", type=str, help="data file")
parser.add_argument("-c", "--coverage_files", nargs="+", help="List of coverage file names")
parser.add_argument("-v", "--stat_files", nargs="+", help="List of stat file names")
parser.add_argument("-m", "--meta_file", type=str, help="metadata file")
parser.add_argument("-r","--real_timecourse", type=func.str2bool)
parser.add_argument("-z","--min_read_count", type=int)
parser.add_argument("-l","--lineage_min_threshold", type=float)
parser.add_argument("-o","--output_file", type=str)
args=parser.parse_args()

# import simulated data + metadata
## read out metadata
meta_df = pd.read_csv(args.meta_file, sep='\t')
meta_df["timepoint"] = meta_df["timepoint"].astype(int)
meta_df["sample_date"] = pd.to_datetime(meta_df["sample_date"], format='%Y-%m-%d')

# keep only data rows from samples list
meta_df = meta_df[meta_df["sample"].isin(args.samples)]

## read out data file
data_df = pd.read_csv(args.data_file, sep=",")
data_df["sample_date"] = pd.to_datetime(data_df["sample_date"], format='%Y-%m-%d')

# keep only data that corresponds with metadata rows. remove unnamed columns
data_df = data_df[data_df["sample_date"].isin(meta_df["sample_date"].to_list())]
data_df = data_df.loc[:, ~data_df.columns.str.contains(r"^Unnamed")]

# create merge dataframe
seq_df = pd.merge(meta_df, data_df, on = "sample_date", how = "outer").sort_values(by=["timepoint"])

# change all lineage columns to float
for col in seq_df.columns:
    if col not in ["sample", "sample_date", "timepoint"]:
        seq_df[col] = seq_df[col].astype(float)

seq_df.replace(0.00, np.nan, inplace=True)

# filter dataframe
seq_df= func.filter_dataframe(seq_df, args.lineage_min_threshold, "simulation")

# # get sequencing quality parameter
seq_df["coverage_avg"] = np.nan
seq_df["coverage_sd"] = np.nan
seq_df["uniformity_wg_per"] = np.nan
seq_df["MAPQ_avg"] = np.nan

for index, row in seq_df.iterrows():
    sample_name = row["sample"]
    i_file = ""

    # get coverage of sample
    for cov_file in args.coverage_files:
        if cov_file.split("/")[-2] == sample_name:
            df = pd.read_csv(cov_file, sep="\t")

            # overall coverage
            seq_df.loc[index, "coverage_avg"] = df[sample_name].mean()
            seq_df.loc[index, "coverage_sd"] = df[sample_name].std()

            # uniformity overall - percentage of bases in all targeted regions (or whole‑genome) 
            # that is covered by at least X%.
            qc_passed = 0

            for i, r in df.iterrows():
                if int(r[sample_name]) >= int(args.min_read_count):
                    qc_passed += 1
            
            seq_df.loc[index, "uniformity_wg_per"] = (qc_passed / len(df)) * 100

            # get average quality score of mapped reads
            for stat_file in args.stat_files:
                if stat_file.split("/")[-2] == sample_name:
                    with open(stat_file) as f:
                        first_line = f.readline()
                        
                    seq_df.loc[index, "MAPQ_avg"] = first_line.split("=")[1].strip()
                    break
            break

# write output file
seq_df.to_csv(args.output_file, index=False)