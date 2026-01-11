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
parser.add_argument("-a", "--abundances_files", nargs="+", help="List of abundance file names")
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

## read out sim data
sim_df = pd.DataFrame({"timepoint": pd.Series(dtype=int)})
variant_list = []

for file in args.abundances_files:
    df = pd.read_csv(file, sep='\t', header=None, names=["lineage", "abundance"])
    
    # check if abundance is 0-100 or 0-1
    if df['abundance'].sum() == 100.0: 
        df['abundance'] = df['abundance'] / 100
        
    df['lineage'] = df['lineage'].str.replace('-cons','')
    variant_list.extend(list(df['lineage']))
    df = df.T
    df.columns = df.iloc[0]
    df = df[1:]

    df["timepoint"] = file[file.rfind("-")+1:file.rfind(".")]
    df["timepoint"] = df["timepoint"].astype(int)
    
    sim_df = pd.concat([sim_df, df], ignore_index=True).sort_values(by=["timepoint"])

# change all columns to float except timepoint
for col in sim_df.columns:
    if col != "timepoint":
        sim_df[col] = sim_df[col].astype(float)
        if not args.real_timecourse:
            sim_df[col] = sim_df[col] / 100

sim_df.replace(0.00, np.nan, inplace=True)
sim_df = sim_df.merge(meta_df, on ="timepoint")

# # get sequencing quality parameter
sim_df["coverage_avg"] = np.nan
sim_df["coverage_sd"] = np.nan
sim_df["uniformity_wg_per"] = np.nan
sim_df["MAPQ_avg"] = np.nan

for index, row in sim_df.iterrows():
    sample_name = row["sample"]
    i_file = ""

    # get coverage of sample
    for cov_file in args.coverage_files:
        if cov_file.split("/")[-2] == sample_name:
            df = pd.read_csv(cov_file, sep="\t")

            # overall coverage
            sim_df.loc[index, "coverage_avg"] = df[sample_name].mean()
            sim_df.loc[index, "coverage_sd"] = df[sample_name].std()

            # uniformity overall - percentage of bases in all targeted regions (or whole‑genome) 
            # that is covered by at least X%.
            qc_passed = 0

            for i, r in df.iterrows():
                if int(r[sample_name]) >= int(args.min_read_count):
                    qc_passed += 1
            
            sim_df.loc[index, "uniformity_wg_per"] = (qc_passed / len(df)) * 100

            # get average quality score of mapped reads
            for stat_file in args.stat_files:
                if stat_file.split("/")[-2] == sample_name:
                    with open(stat_file) as f:
                        first_line = f.readline()
                        
                    sim_df.loc[index, "MAPQ_avg"] = first_line.split("=")[1].strip()
                    break
            break

# filter dataframe
sim_df= func.filter_dataframe(sim_df, args.lineage_min_threshold, "simulation")

# write output file
sim_df.to_csv(args.output_file, index=False)