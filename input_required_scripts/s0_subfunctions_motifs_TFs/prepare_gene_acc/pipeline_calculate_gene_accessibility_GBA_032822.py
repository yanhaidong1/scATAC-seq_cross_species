#!/usr/bin/env python

##updating 032822 only consider the GBM not consider the cicero
##this script will generate a pipeline to calculate the gene accessibility

import re
import glob
import sys
import subprocess
import os


input_all_required_pipeline_scripts_dir = sys.argv[1]

input_meta_fl = sys.argv[2]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/08_training_classifiers_031822/raw_data/Arab_root/Arabidopsis_root_scATAC_atlas.annotation.txt

tissue_name = sys.argv[3] ##All

#########
##step 02 calculate the accessibility of the normGBA
input_gene_sparse_fl = sys.argv[4]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/08_training_classifiers_031822/raw_data/Arab_root/At_root.all.genes.sparse

##/scratch/hy17471/rice_mPing_scATAC_seq_050521/07_1_prepare_gene_acc_files_060421/output_dir/step03_prepare_gene_sparse_dir/02_TPO_norm_dir/opt_gene_sparse_TPO_norm.sparse

input_gene_bed_fl = sys.argv[5]
##/scratch/hy17471/rice_mPing_scATAC_seq_050521/07_1_prepare_gene_acc_files_060421/output_dir/01_generate_gene_bed_fl_dir/opt_riceMSUr7_genes_500bpTSS_sorted.bed

input_output_dir = sys.argv[6]




##define a function for the gene body sparse
def step02_add_gene_body_acc (input_all_required_pipeline_scripts_dir,input_gene_sparse_fl,input_meta_fl,input_gene_bed_fl,tissue_name,output_dir):

    call_normGBA_script = input_all_required_pipeline_scripts_dir + '/02_call_accessiblity_normGBA.R'

    cmd = 'Rscript ' + call_normGBA_script + \
          ' ' + input_gene_sparse_fl + \
          ' ' + input_meta_fl + \
          ' ' + input_gene_bed_fl + \
          ' ' + tissue_name + \
          ' ' + output_dir
    subprocess.call(cmd,shell=True)


#########
##step 02
store_call_accessiblity_normGBA_dir = input_output_dir + '/02_call_accessiblity_normGBA_dir'
if not os.path.exists(store_call_accessiblity_normGBA_dir):
    os.makedirs(store_call_accessiblity_normGBA_dir)

step02_add_gene_body_acc (input_all_required_pipeline_scripts_dir,input_gene_sparse_fl,input_meta_fl,input_gene_bed_fl,tissue_name,store_call_accessiblity_normGBA_dir)




