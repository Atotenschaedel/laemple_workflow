# -*- coding: utf-8 -*-
"""
Created on 29-06-2023
This script generates meta data file for vaquero tool based on parameter from config file
@author: aschedl
"""
# import libraries
import argparse
import pandas as pd

# add argument
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str)
parser.add_argument("-o","--output", type=str)
parser.add_argument("--run", type=str)
parser.add_argument("--start_date", type=str)
parser.add_argument("--locationID", type=str)
parser.add_argument("--locationName", type=str)
parser.add_argument("--N_in_consensus", type=str)
parser.add_argument("--adress_town", type=str)
parser.add_argument("--connected_people", type=str)
parser.add_argument("--latitude", type=str)
parser.add_argument("--longitude", type=str)
parser.add_argument("--include_in_report", type=str)
parser.add_argument("--report_category", type=str)
parser.add_argument("--status", type=str)
parser.add_argument("--suffix", type=str)

args=parser.parse_args()

# create dataframe
df = pd.read_csv(args.input, sep="\t")

df["BSF_run"] = args.run
df["BSF_start_date"] = args.start_date
df["LocationID"] = args.locationID
df["LocationName"] = args.locationName
df["N_in_Consensus"] = args.N_in_consensus
df["adress_town"] = args.adress_town
df["connected_people"] = args.connected_people
df["dcpLatitude"] = args.latitude
df["dcpLongitude"] = args.longitude
df["include_in_report"] = args.include_in_report
df["report_category"] = args.report_category
df["status"] = args.status
df["additional_information"] = ""
df["BSF_sample_name"] = df["sample"] + args.suffix

df.rename(columns={"sample":"RNA_ID_int"}, inplace=True)

df.to_csv(args.output, index=False, sep="\t")