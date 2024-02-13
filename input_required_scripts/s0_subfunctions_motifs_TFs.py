#!/usr/bin/env python


##Here we will check the TF accessiblity

import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats
from Bio import SeqIO
from statistics import mean



########
##step01 prepare the gene accessiblity
def subfunction_prepare_gene_accessibility (input_gene_cell_raw_sparse_fl,input_meta_fl,input_required_scripts_dir,spe_prefix,
                                            opt_dir):

    call_normGBA_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/prepare_gene_acc/02_call_accessiblity_normGBA_notrequire_genebed.R'

    cmd = 'Rscript ' + call_normGBA_script + \
          ' ' + input_gene_cell_raw_sparse_fl + \
          ' ' + input_meta_fl + \
          ' ' + spe_prefix + \
          ' ' + opt_dir
    print(cmd)
    subprocess.call(cmd, shell=True)


def subfunction_prepare_smooth_gene_accessibility (input_required_scripts_dir,input_meta_fl,input_gene_norm_sparse_fl,input_svd_fl,
                                                    input_core_num,s0_s5_sub_target_cluster,
                                                   opt_dir):


    ##the first round will run all the markers at the first time
    ##if we obtain file from the first round we will consider use the opt to generate markers as we want
    gene_smooth_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/smooth/smooth_genes_acc_112723.R'
    function_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/smooth/subfunction_smooth_genes_acc_112723.R'

    cmd = 'Rscript ' + gene_smooth_script + \
          ' ' + input_meta_fl + \
          ' ' + input_gene_norm_sparse_fl + \
          ' ' + input_svd_fl + \
          ' ' + input_core_num + \
          ' ' + s0_s5_sub_target_cluster + \
          ' ' + opt_dir + \
          ' ' + function_script
    print(cmd)
    subprocess.call(cmd, shell=True)


########
##step02 we will smooth motif
def subfunction_prepare_smooth_motif_deviation (input_required_scripts_dir,ipt_meta_fl,ipt_motif_deviation_fl,
                                                ipt_svd_fl,spe_prefix,s0_s5_sub_target_cluster,opt_dir):


    motif_smooth_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/smooth/smooth_motif_deviation_112723.R'

    print('- open conduct smooth for the motif deviation score')
    cmd = 'Rscript ' + motif_smooth_script + \
          ' ' + ipt_meta_fl + \
          ' ' + ipt_motif_deviation_fl + \
          ' ' + ipt_svd_fl + \
          ' ' + spe_prefix + \
          ' ' + s0_s5_sub_target_cluster + \
          ' ' + opt_dir
    print(cmd)
    subprocess.call(cmd, shell=True)

    print('- open transfer smooth sparse')
    motif_smooth_transfer_format_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/smooth/smooth_motif_deviation_transfer_format_112723.R'
    smoothed_motifs_fl = opt_dir + '/' + spe_prefix + '.smoothed_motifs.txt'
    cmd = 'Rscript ' + motif_smooth_transfer_format_script + \
          ' ' + smoothed_motifs_fl + \
          ' ' + opt_dir
    print(cmd)
    subprocess.call(cmd, shell=True)


########
##step03 we will correlate the TF and motif
def subfunction_correlate_TF_motif (input_required_scripts_dir,ipt_motif_smooth_sparse_fl,ipt_meta_fl,
                                    ipt_smooth_TF_mtx_rds_fl,ipt_correspond_At_gene_to_spe_fl,input_At_TF_to_motif_fl,input_Atgene_TFid_fl,
                                    s0_s5_sub_target_cluster,opt_dir):

    correlate_TF_motif_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/correlate_motifs_TFs/correlation_motif_TFs.py'
    subfunction_correlation_R_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/correlate_motifs_TFs/subfunction_correlation.R'

    ###########################
    ##for the motif avg cluster
    ipt_smooth_sparse_fl = ipt_motif_smooth_sparse_fl
    ipt_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/correlate_motifs_TFs/subfunction_avg_smooth_motif_dev.R'

    cmd = 'Rscript ' + ipt_script + \
          ' ' + ipt_smooth_sparse_fl + \
          ' ' + ipt_meta_fl + \
          ' ' + s0_s5_sub_target_cluster + \
          ' ' + opt_dir
    print(cmd)
    subprocess.call(cmd, shell=True)

    norm_motif_avg_cluster_fl = opt_dir + '/opt_motif_average_smooth_accessible_clusters.txt'


    ########################
    ##for the TF avg cluster
    ipt_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/correlate_motifs_TFs/subfunction_avg_smooth_TF_dev.R'

    cmd = 'Rscript ' + ipt_script + \
          ' ' + ipt_smooth_TF_mtx_rds_fl + \
          ' ' + ipt_meta_fl + \
          ' ' + s0_s5_sub_target_cluster + \
          ' ' + opt_dir
    print(cmd)
    subprocess.call(cmd, shell=True)

    input_norm_downsampling_gene_act_final_fl = opt_dir + '/opt_TF_average_smooth_accessible_clusters.txt'

    #########################
    ##conduct the correlation
    ##step02_s4_3_sub_open_compare_corrRate will be set as no
    ##we only care about the max correlation
    cmd = 'python ' + correlate_TF_motif_script + \
          ' ' + ipt_correspond_At_gene_to_spe_fl + \
          ' ' + input_Atgene_TFid_fl + \
          ' ' + input_At_TF_to_motif_fl + \
          ' ' + input_norm_downsampling_gene_act_final_fl + \
          ' ' + norm_motif_avg_cluster_fl + \
          ' ' + subfunction_correlation_R_script + \
          ' ' + 'no' + \
          ' ' + opt_dir
    print(cmd)
    subprocess.call(cmd, shell=True)



















