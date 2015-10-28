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
    if data[0] == 'symbol':
        continue
    symbol = data[0]
    rpm = data[1]
    ref_dict[symbol] = rpm

input_file = open(PATH_INPUT_FILE, 'r')
output_file = open(PATH_OUTPUT_FILE, 'w')

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == "gr_id":
        print('gr_id','symbol','refid_representative','seq_length','rpm','rpkm',sep="\t",end="\n",file=output_file)
        continue
    gr_id = data[0]
    symbol = data[1]
    refid_representative = data[2]
    seq_length = data[3]
    rpm = ref_dict[symbol]
    rpkm = 0
    if seq_length == 'NA':
        pass
    else:
        rpkm = float(rpm)/float(seq_length)*1000
    print(gr_id,symbol,refid_representative,seq_length,rpm,rpkm, sep="\t",end="\n",file=output_file)
