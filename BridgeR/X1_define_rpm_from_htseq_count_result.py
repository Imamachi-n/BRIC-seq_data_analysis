#!/usr/bin/env python3

import re
import sys

PATH_REF_FILE = sys.argv[1]
PATH_INPUT_FILE = sys.argv[1]
PATH_OUTPUT_FILE = sys.argv[2]

#ref_file = open('/mnt/hgfs/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_siStealth_siPUM1_htseq_count/htseq_count_result_siPUM1_0h.txt','r')
ref_file = open(PATH_REF_FILE,'r')

total_count = 0

for line in ref_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == '' or re.match('^__',data[0]):
        continue
    symbol = data[0]
    read_count = int(data[1])
    total_count += read_count

ref_file.close()

#input_file = open('/mnt/hgfs/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_siStealth_siPUM1_htseq_count/htseq_count_result_siPUM1_0h.txt','r')
#output_file = open('/mnt/hgfs/github/BRIC-seq_data_analysis/BridgeR/data/BridgeR_siStealth_siPUM1_htseq_count/htseq_count_result_siPUM1_0h_rpm.txt','w')
input_file = open(PATH_INPUT_FILE,'r')
output_file = open(PATH_OUTPUT_FILE,'w')

print('symbol','rpm',sep="\t",end="\n",file=output_file)

for line in input_file:
    line = line.rstrip()
    data = line.split("\t")
    if data[0] == '' or re.match('^__',data[0]):
        continue
    symbol = data[0]
    rpm = int(data[1])/total_count*1000000
    print(symbol,rpm,sep="\t",end="\n",file=output_file)

input_file.close()
output_file.close()
