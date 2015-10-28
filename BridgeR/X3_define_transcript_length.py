#!/usr/bin/env python3

import re
import sys

PATH_INPUT_FILE = sys.argv[1]
PATH_OUTPUT_FILE = sys.argv[2]

fasta_file = open(PATH_INPUT_FILE,'r')
output_file = open(PATH_OUTPUT_FILE, 'w')
fasta_dict = {}
checker = []
name = ''
tmp2 = ''
for line in fasta_file:
    line = line.rstrip()
    if line.startswith('#'):
        continue
    if line.startswith('>'):
        name = line[1:].strip()
        if not name in checker:
            checker.append(name)
        else:
            print ('ERROR: The same name exists =>' + name)
            sys.exit(1)
            continue
        fasta_dict[name] = 0
    else:
        fasta_dict[name] += len(line)
fasta_file.close()

for x in fasta_dict.keys():
    seq_length = fasta_dict[x]
    print(x,seq_length, sep="\t", end="\n", file=output_file)
