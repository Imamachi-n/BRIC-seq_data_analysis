#!/usr/bin/env python3
#Title: A2_merge_bowtie24Luc2_and_cuffnorm4rpkm_infor
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-21

import re
output_file = open("BRIC-seq_PUM1_study_Luc2_Rel_RPKM.txt",'w')

def bowtie2_input(input_file_infor):
    #bowtie2 - result
    #4417 (0.44%) aligned exactly 1 time
    regex = r'\s+(?P<reads>[0-9]+) \([0-9]+\.[0-9]+%\) aligned exactly 1 time'
    reads_list = []

    for input_files in input_file_infor:
        input_file = open(input_files,'r')
        for line in input_file:
            line = line.rstrip()
            if re.search('aligned exactly 1 time',line):
                result = re.match(regex, line)
                reads = int(result.group('reads'))
                if not reads_list:
                    reads_list = [reads]
                else:
                    reads_list.append(reads)
    return reads_list

def htseq_count_total_reads(input_count_file_infor):
    total_reads_list = []
    for input_files in input_count_file_infor:
        input_file = open(input_files,'r')
        total_count = 0
        for line in input_file:
            line = line.rstrip()
            data = line.split("\t")
            if data[0] == '' or re.match('^__',data[0]):
                continue
            count_reads = data[1]
            total_count += int(count_reads)
        if not total_reads_list:
            total_reads_list = [total_count]
        else:
            total_reads_list.append(total_count)
    return total_reads_list

#################################################
print('sample','T0_1','T1_1','T2_1','T4_1','T8_1','T12_1',sep="\t",end="\n",file=output_file)


########################
input_file_infor = ["unmapped_bowtie2_infor_BRIC-seq_siStealth_0h.txt",
                    "unmapped_bowtie2_infor_BRIC-seq_siStealth_1h.txt",
                    "unmapped_bowtie2_infor_BRIC-seq_siStealth_2h.txt",
                    "unmapped_bowtie2_infor_BRIC-seq_siStealth_4h.txt",
                    "unmapped_bowtie2_infor_BRIC-seq_siStealth_8h.txt",
                    "unmapped_bowtie2_infor_BRIC-seq_siStealth_12h.txt"]

input_count_file_infor = ["htseq_count_result_siStealth_0h.txt",
                          "htseq_count_result_siStealth_1h.txt",
                          "htseq_count_result_siStealth_2h.txt",
                          "htseq_count_result_siStealth_4h.txt",
                          "htseq_count_result_siStealth_8h.txt",
                          "htseq_count_result_siStealth_12h.txt"]

reads_list = bowtie2_input(input_file_infor)
reads_list = list(map(float,reads_list))
total_reads_list = htseq_count_total_reads(input_count_file_infor)
total_reads_list = list(map(float,total_reads_list))

rel_rpkm = [reads_list[x]/total_reads_list[x]*1000000 for x in range(0,len(reads_list))]
rel_rpkm = [rel_rpkm[x]/rel_rpkm[0] for x in range(0,len(rel_rpkm))]
rel_rpkm = list(map(str,rel_rpkm))
print('siStealth',"\t".join(rel_rpkm),sep="\t",end="\n",file=output_file)

########################

input_file_infor = ["unmapped_bowtie2_infor_BRIC-seq_siPUM1_0h.txt",
                    "unmapped_bowtie2_infor_BRIC-seq_siPUM1_1h.txt",
                    "unmapped_bowtie2_infor_BRIC-seq_siPUM1_2h.txt",
                    "unmapped_bowtie2_infor_BRIC-seq_siPUM1_4h.txt",
                    "unmapped_bowtie2_infor_BRIC-seq_siPUM1_8h.txt",
                    "unmapped_bowtie2_infor_BRIC-seq_siPUM1_12h.txt"]
input_count_file_infor = ["htseq_count_result_siPUM1_0h.txt",
                          "htseq_count_result_siPUM1_1h.txt",
                          "htseq_count_result_siPUM1_2h.txt",
                          "htseq_count_result_siPUM1_4h.txt",
                          "htseq_count_result_siPUM1_8h.txt",
                          "htseq_count_result_siPUM1_12h.txt"]

reads_list = bowtie2_input(input_file_infor)
reads_list = list(map(float,reads_list))

total_reads_list = htseq_count_total_reads(input_count_file_infor)
total_reads_list = list(map(float,total_reads_list))

rel_rpkm = [reads_list[x]/total_reads_list[x]*1000000 for x in range(0,len(reads_list))]
rel_rpkm = [rel_rpkm[x]/rel_rpkm[0] for x in range(0,len(rel_rpkm))]
rel_rpkm = list(map(str,rel_rpkm))
print('siPUM1',"\t".join(rel_rpkm),sep="\t",end="\n",file=output_file)

