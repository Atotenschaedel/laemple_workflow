#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 11:18:23 2021

@author: lendler
"""

from pysam import VariantFile
import sys
import argparse

parser = argparse.ArgumentParser(description="""
takes a vcf file as created by freebayes and calculates and adds the format/AF field
if the input is stdin, just use "-" for input (default)
writes to std out
""")
parser.add_argument("-i", dest="vcf_file", help="vcf input file (can be gz), for stdin use -", default="-")
args = parser.parse_args()
#vcf_file = "/Users/lendler/HPC_lab_bergthaler/2020_SARS-CoV-2_Evolution/run_BSF_0935_H5K2CDRXY/SarsVirSeq2/results_pipeline/CoV_5963_S80668/CoV_5963_S80668_samp.lofreq_filtered_phased.vcf.gz"
#vcf_file = "/Users/lendler/Data/Test_virseq/CoV_13511_S93769/CoV_13511_S93769_samp.lofreq_filtered_phased.vcf.gz"

vcf_in = VariantFile(args.vcf_file)
# check if AF missing in fields
if 'AF' in vcf_in.header.formats.keys():
    sys.exit("Already an AF field in formats!")
elif len(set(["DP", "AD"]).intersection(vcf_in.header.formats.keys())) < 2:
    sys.exit("Missing AD and/or DP fields in formats!")

# add afs
vcf_in.header.formats.add("AF", "A", "Float", "Allele frequency derived from reads")

# create output file with new header
vcf_out = VariantFile("-", mode='w', header=vcf_in.header)
#print(vcf_in.header)

for rec in vcf_in:
    for sample in rec.samples:
        rec.samples[sample]['AF'] = [float(x) / rec.samples[sample]['DP']
                                     for x in rec.samples[sample]['AD'][1:len(rec.samples[sample]['AD'])]]
    #print(rec)
    vcf_out.write(rec)

vcf_out.close()
