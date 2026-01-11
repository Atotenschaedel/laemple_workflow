# -*- coding: utf-8 -*-
"""
Created on Sun May 21 10:30:22 2023
This script is doing ....
@author: aschedl
"""
from Bio import SeqIO
import os

w_dir = "experiments/Ex01_WideHighQual/data"

f_list=[]

for file in os.listdir(w_dir):
    f = os.path.join(w_dir, file)
    f_list.append(f)

def check(test_seq):
    # allow IUPAC nucleotide code expect gap (-/.)
    allowed=set("ACGTURYSWKMBDHVN")
    return set(test_seq) <= allowed


print("we are starting")

for index, fastq in enumerate(f_list):
    file = open(fastq, 'r')
    print(fastq)
    problem_lines = []

    count = 0

    while True:
        count += 1
        line = file.readline().strip()
        if count % 2 == 0 and not count % 4 == 0: # is sequence
            if not check(line):
            #if "-" in line:
                #print("Line{}: Sequence: {}".format(count, line))
                indices = [index for index in range(len(line)) if line.startswith("-", index)]
                #print("found special character at: ", indices)
                problem_lines.append([count, indices])
                
                #print("corrected string:")
                #print("".join([char for idx, char in enumerate(line) if idx not in indices]))
        
#           elif count % 4 == 0:
#              print("Line{}: Quality: {}".format(count, line.strip()))

        if not line:
            break
    
    file.close()

    if len(problem_lines) > 0:
        print("we have problems in: ", file)