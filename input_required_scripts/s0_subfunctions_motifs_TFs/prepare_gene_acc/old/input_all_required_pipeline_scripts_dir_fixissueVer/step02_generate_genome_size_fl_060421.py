#!/usr/bin/env python

##this script we will combine the root hair and non root hair data
import re
import glob
import os
import sys
import subprocess
import random
from Bio import SeqIO

##we will generate the genome size fl
input_genome_fl = sys.argv[1]
input_output_dir = sys.argv[2]

def generate_genome_size (input_genome_fl,input_output_dir):

    store_final_seq_dic = {}
    for seq_recorde in SeqIO.parse(input_genome_fl,'fasta'):
        store_final_seq_dic[seq_recorde.id] = str(len(str(seq_recorde.seq)))

    store_final_line_list = []
    for eachid in store_final_seq_dic:
        final_line = eachid + '\t' + store_final_seq_dic[eachid]
        store_final_line_list.append(final_line)

    with open (input_output_dir + '/opt_genome_size.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

generate_genome_size (input_genome_fl,input_output_dir)
