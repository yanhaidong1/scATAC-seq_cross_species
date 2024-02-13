#!/usr/bin/env python

##updating 060421 consider the rice case
##updating 012921 build the nonbinary gene sparse fl
##this script will help us to generate the gene sparse fl based on the temp fl

import re
import glob
import sys
import subprocess
import os

input_temp_unique_bed_fl = sys.argv[1] ##00_generate_unique_TSS_BC_counts_110420/temp_unique_bed_dir
input_gene_bed_fl = sys.argv[2]##resources/add_00_prepare_gene_bed_fl_120220/opt_Wm82a4v1_genes_500bpTSS_sorted.bed
input_genome_fai_fl = sys.argv[3] ##resources/genome.fa.fai.sorted
input_output_dir = sys.argv[4] ##resources/add_00_prepare_gene_sparse_fl_120220/output_dir

input_fastSparsetn5_pl = sys.argv[5]
##scripts/fastSparse.tn5.pl ##change to non-binary counts fastSparse.nonbinary.peak.pl
##scripts/fastSparse.nonbinary.peak.pl

def make_sparse (input_temp_unique_bed_fl,input_gene_bed_fl,input_genome_fai_fl,input_output_dir,input_fastSparsetn5_pl):

    opt_gene_sparse_bed_dir = input_output_dir + '/opt_gene_sparse_bed_dir'
    if not os.path.exists(opt_gene_sparse_bed_dir):
        os.makedirs(opt_gene_sparse_bed_dir)

    cmd = 'bedtools intersect -a ' +  input_temp_unique_bed_fl + \
          ' -b ' + input_gene_bed_fl + ' -wa -wb' + \
          ' -g ' + input_genome_fai_fl + \
          ' -sorted | perl ' + input_fastSparsetn5_pl + \
          ' - > ' + opt_gene_sparse_bed_dir + '/opt_gene_chnm.sparse'
    subprocess.call(cmd,shell=True)

def change_name_gene (opt_combine_gene_fl,input_output_dir):

    ##change the location name to the gene
    store_loc_gene_dic = {}
    with open (input_gene_bed_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            gene_loc = col[0] + '_' + col[1] + '_' + col[2]
            gene_name = col[3]
            store_loc_gene_dic[gene_loc] = gene_name

    ##the following script will allow large mem we will directly add information for eachline wehen we read the sam file
    ipt_gene_sparse_fl = open(opt_combine_gene_fl, 'r')
    opt_chnm_fl = open(input_output_dir + '/opt_gene_sparse_bed_dir/opt_gene_chnm.sparse', 'w')  # w when u wanna write sth on the file

    store_wrong_loc_dic = {}
    for eachline in ipt_gene_sparse_fl:
        eachline = eachline.strip('\n')
        col = eachline.strip().split()
        if col[0] in store_loc_gene_dic:
            final_line = store_loc_gene_dic[col[0]] + '\t' + col[1] + '\t' + col[2]
            opt_chnm_fl.write(final_line + '\n')
        else:
            store_wrong_loc_dic[col[0]] = 1

    ipt_gene_sparse_fl.close()
    opt_chnm_fl.close()

    with open (input_output_dir + '/opt_gene_sparse_bed_dir/opt_wrong_loc.txt','w+') as opt:
        for eachline in store_wrong_loc_dic:
            opt.write(eachline + '\n')


make_sparse (input_temp_unique_bed_fl,input_gene_bed_fl,input_genome_fai_fl,input_output_dir,input_fastSparsetn5_pl)
#change_name_gene (input_output_dir + '/opt_gene_sparse_bed_dir/opt_gene.sparse',input_output_dir)
