#!/usr/bin/env python3

import re
import sys

PATH_REF_FILE = sys.argv[1]
PATH_INPUT_FILE = sys.argv[2]
PATH_OUTPUT_FILE = sys.argv[3]

ref_file = open(PATH_REF_FILE,'r')

ref_dict = {}

for line in ref_file:
    line = line.rstrip()
    data = line.split("\t")
    refid = data[0]
    seq_length = data[1]
    ref_dict[refid] = seq_length

input_file = open(PATH_INPUT_FILE, 'r')
output_file = open(PATH_OUTPUT_FILE, 'w')

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == "gr_id":
        print('gr_id','symbol','refid_representative','seq_length',sep="\t",end="\n",file=output_file)
        continue
    seq_length = 'NA'
    if data[2] in ref_dict:
        seq_length = ref_dict[data[2]]
    else:
        seq_length = "NA"
    print(data[0],data[1],data[2],seq_length, sep="\t",end="\n",file=output_file)
