# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 10:30:22 2015
This script uses pysam and bedtools to downsample a bam file to an approximate coverage threshold.
All reads are shuffled and added until their median coverage for the region covered by each read surpasses the threshold.
This can lead to quite some variation around the threshold, although not in the order of magnitudes. 
The algorithm is quite inefficient and needs ridiculous amounts of memory ( around 8.5 GB per 32 Mio reads), as it store all read names in a hash and a list.
It gets the intervals for reads mapping from bam2bed from bedtools (should be in path, else give path to it)
@author: lukasendler
Update 23-05-2023, added additional parameter "conda-env", if true python script can be executed inside a shell with all bedtools already installed 
@author: aschedl
"""

def get_read_intervals(bam_file,bedtools,verb, conda):
    """ 
    get a bam file name, the path to bam2bed and a verbose flag
    pipes the ouput of bedtools bamtobed and creates a dictionary of the reads with covered regions 
    returns dict: chrom names -> idx
    dict: readname -> [idx, start, end]    
    """
    read_dict=defaultdict(list)
    chroms=defaultdict()
    mateid=re.compile("[/_][12]$")
    if conda:
        p=subprocess.Popen(f"bedtools bamtobed -i {bam_file}", shell=True,stdout=subprocess.PIPE,bufsize=1024)
    else:
        p=subprocess.Popen([ bedtools , "bamtobed" , "-i" ,bam_file ],stdout=subprocess.PIPE,bufsize=1024)
    lc=0
    for line in p.stdout:
        if verb:
            lc += 1
            if lc % 1000000 == 0:
                sys.stderr.write("loaded read coverage intervals: " + str(lc) +"\n")
                sys.stderr.flush()
        # get read info, need to decode as otherwise bit strings
        entries=line.decode('utf8').split()[0:4]
        # if new chromosome, add to chroms dict
        if not entries[0] in chroms.keys():
            chroms[entries[0]] = len(chroms.keys())
        # append covered region to entries
        read_dict[re.sub(mateid,"",entries[3])].append([ chroms[entries[0]], int(entries[1]),int(entries[2]) ])            
    return(chroms, read_dict)

def check_cov(cov_arrays,covs,max_cov,start_only):
    """
    gets an list of lists with chrom idxs and the start & ends of coverages, the coverage arrays and the max coverage
    checks whether coverage medians are below the max threshold in the covered regions
    returns True if yes
    also add values to the covered streches in covs
    if start_only, only check and add coverage to start of read
    """
    # get median values of covered regions 
    if start_only:
        ret = True
        for x in cov_arrays:
            if covs[x[0]][x[1]] > max_cov:
                ret = False
            covs[x[0]][x[1]] += 1
        return ret
    min_med_cov=np.min([ np.median(covs[x[0]][x[1]:x[2]]) for x in cov_arrays])
    if min_med_cov < max_cov:
        for x in cov_arrays:
            # add 1 to covered regions
            covs[x[0]][x[1]:x[2]] += 1
        return True
    else:
        return False
    
import pysam
#import gzip
import numpy as np
import sys, os, time
import subprocess
import re
#from Bio import trie
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="""
gets a bam file name and a coverage threshold and tries to a give a randomly sampled bam file with a max. coverage around the threshold.
The 
""") 
parser.add_argument("-b", dest="bam_file", help="bam file", required=True)
parser.add_argument("-c", dest="max_cov", help="maximal coverage threshold (def.: 10000)", default=10000)
parser.add_argument("--bedtools", dest="bedtools", help="executable bedtools with path (def.: \"/usr/local/bin/bedtools\")", default="/usr/local/bin/bedtools")
parser.add_argument("--maxlen", dest="maxlen", help="maximal length of chromosomes  (def.: 100000)", default=100000)
parser.add_argument("--start", dest="start_only", action="store_true", help="only check coverage at 5' of reads, similar to GATK, for amplicons (def.: FALSE)", default=False)
parser.add_argument("-v","--verb", dest="verb", action="store_true", help="print verbose state on stderr", default=False)
parser.add_argument("--conda-env", dest="conda", action="store_true", help="execute script inside conda environment with bedtools and installed", default=False)

args = parser.parse_args()
bam_file = vars(args)['bam_file']
max_cov = int(vars(args)['max_cov'])
maxlen = int(vars(args)['maxlen'])
bedtools = vars(args)['bedtools']
verb = vars(args)['verb']
start_only = vars(args)['start_only']
conda=vars(args)['conda']

verb=True
if verb:
    sys.stderr.write("at " + time.strftime("%c") + " starting with max_cov:" + str(max_cov) + " bamfile:" + bam_file + "\n")
# start bam2bedpipe and fill up readname and interval dict
(chroms, read_dict)=get_read_intervals(bam_file,bedtools,verb,conda)
#create coverage arrays
covs=[ np.zeros(maxlen,dtype=int) for x in chroms.keys()]
# get readname list and shuffle it
read_names=[ x for x in read_dict.keys()]
shufs=np.arange(0,len(read_names))
np.random.shuffle(shufs)
# check which reads to keep
keeps=[]
if verb:
    sys.stderr.write("at " + time.strftime("%c") + " need to check coverages for " + str(len(read_names)) + " reads\n")
lc=0
for i in shufs:
    if verb:
        lc += 1
        if lc % 1000000 == 0:
            sys.stderr.write("checked coverage of "+ str(lc) + " readpairs\n")
            sys.stderr.flush()
    if check_cov(read_dict[read_names[i]],covs,max_cov,start_only):
        keeps.append(read_names[i])

# remove read_dict, shufs and read_names
del read_dict
del read_names
del shufs
# transform to read set (for fast lookup)
keeps=set(keeps)

# read bam file read by read
samfile=pysam.Samfile(bam_file,"rb")
outf = bam_file
# remove .bam at end
if outf[-4:] == ".bam":
    outf = outf[:-4]
outsam=pysam.Samfile(outf+"_max_cov_"+str(max_cov)+".bam","wb",template=samfile)
## split BAM file
if verb:
    sys.stderr.write("at " + time.strftime("%c") + "writing reduced coverage file\n")
lc=0
for l in samfile:
    if verb:
        lc += 1
        if lc % 1000000 == 0:
            sys.stderr.write("went thorough "+ str(lc) + " reads\n")
            sys.stderr.flush()
    if l.qname in keeps:
        # if in set, write
        outsam.write(l)
outsam.close()
if verb:
    sys.stderr.write("at " + time.strftime("%c") + "finished writing\n")
samfile.close() 


 
