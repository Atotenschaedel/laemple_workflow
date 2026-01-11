# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 10:30:22 2023
@author: aschedl
Description: This script takes a tsv file with coverage per positon and saves plot as png
used via command line
"""
# import libraries
import matplotlib.pyplot as plt
import pandas as pd
import argparse

# add arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input_file", type=str)
parser.add_argument("-o","--output_file", type=str)
parser.add_argument("-s","--sample_name", type=str)
args=parser.parse_args()

# read out file + plot
df = pd.read_csv(args.input_file, sep="\t")
df.plot(x="POS", y=args.sample_name)
plt.savefig(args.output_file)
plt.close()