#!/usr/bin/env python

##updating 120423 add a function to plot the family

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
from operator import itemgetter



##plot motifs or TFs
def subfunction_plot_TFs (input_required_scripts_dir,ipt_TF_motif_corr_fl,ipt_gene_impute_acc_rds_fl,GAobj_rds_fl
                          ,opt_dir,s0_s5_sub_lim,s0_s5_sub_plot_topnum_pairs):

    ##opt_corresponding_motif_TF.txt is the ipt_TF_motif_corr_fl
    ##s0_s5_sub_plot_topnum_pairs is the number of
    store_TF_motif_pair_coef_dic = {}
    with open (ipt_TF_motif_corr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            motif_TF = col[0] + '__' + col[1]
            coef = float(col[2])
            store_TF_motif_pair_coef_dic[motif_TF] = coef

    ##select the top for building the candidate markers
    res = dict(sorted(store_TF_motif_pair_coef_dic.items(), key=itemgetter(1), reverse=True)[:int(s0_s5_sub_plot_topnum_pairs)])

    store_target_gene_dic = {}
    for eachpair in res:
        mt = re.match('(.+)__(.+)',eachpair)
        motif = mt.group(1)
        gene = mt.group(2)
        store_target_gene_dic[gene] = 1

    store_final_line_list = []
    first_line = 'name' + '\t' + 'geneID' + '\t' +  'type' + '\t' + 'tissue' + '\t' + 'common'
    store_final_line_list.append(first_line)

    for eachgene in store_target_gene_dic:
        final_line = eachgene + '\t' + eachgene + '\t' + 'TF' + '\t' + 'all' + '\t' + eachgene
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_top_' + s0_s5_sub_plot_topnum_pairs + '_TF_markers.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    plot_TFs_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/plot_TFs_motifs/plot_genes_112723.R'

    cmd = 'Rscript ' + plot_TFs_script + \
          ' ' + ipt_gene_impute_acc_rds_fl + \
          ' ' + GAobj_rds_fl + \
          ' ' + opt_dir + '/opt_top_' + s0_s5_sub_plot_topnum_pairs + '_TF_markers.txt' + \
          ' ' + opt_dir + \
          ' ' + s0_s5_sub_lim
    print(cmd)
    subprocess.call(cmd,shell=True)


##plot motifs or TFs
def subfunction_plot_motifs (input_required_scripts_dir,ipt_TF_motif_corr_fl,ipt_motif_impute_acc_threesparse_fl, ipt_meta_fl,
                          opt_dir,s0_s5_sub_lim,s0_s5_sub_plot_topnum_pairs):

    ##opt_corresponding_motif_TF.txt is the ipt_TF_motif_corr_fl
    ##s0_s5_sub_plot_topnum_pairs is the number of
    store_TF_motif_pair_coef_dic = {}
    with open (ipt_TF_motif_corr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            motif_TF = col[0] + '__' + col[1]
            coef = float(col[2])
            store_TF_motif_pair_coef_dic[motif_TF] = coef

    ##select the top for building the candidate markers
    res = dict(sorted(store_TF_motif_pair_coef_dic.items(), key=itemgetter(1), reverse=True)[:int(s0_s5_sub_plot_topnum_pairs)])

    store_target_motif_dic = {}
    for eachpair in res:
        mt = re.match('(.+)__(.+)',eachpair)
        motif = mt.group(1)
        gene = mt.group(2)
        store_target_motif_dic[motif] = 1

    store_final_line_list = []
    first_line = 'name' + '\t' + 'geneID' + '\t' +  'type' + '\t' + 'tissue' + '\t' + 'common'
    store_final_line_list.append(first_line)

    for eachmotif in store_target_motif_dic:
        final_line = eachmotif + '\t' + eachmotif + '\t' + 'TF' + '\t' + 'all' + '\t' + eachmotif
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_top_' + s0_s5_sub_plot_topnum_pairs + '_motif_markers.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    plot_motifs_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/plot_TFs_motifs/plot_motifs_112723.R'

    cmd = 'Rscript ' + plot_motifs_script + \
          ' ' + ipt_motif_impute_acc_threesparse_fl + \
          ' ' + ipt_meta_fl + \
          ' ' + opt_dir + '/opt_top_' + s0_s5_sub_plot_topnum_pairs + '_motif_markers.txt' + \
          ' ' + opt_dir + \
          ' ' + s0_s5_sub_lim
    print(cmd)
    subprocess.call(cmd,shell=True)

##updating 120423
def subfunction_plot_TFfam_motifs_TFs (input_required_scripts_dir,ipt_TF_motif_corr_fl,ipt_motif_commonnm_fl,
                                       ipt_gene_impute_acc_rds_fl, GAobj_rds_fl,
                                       ipt_motif_impute_acc_threesparse_fl,ipt_meta_fl,
                                       opt_dir,s0_s5_sub_target_fam_str,s0_s5_sub_lim):

    s0_s5_sub_target_fam_list = s0_s5_sub_target_fam_str.split(',')

    for eachfam in s0_s5_sub_target_fam_list:

        store_TFfam_dir = opt_dir + '/' + eachfam
        if not os.path.exists(store_TFfam_dir):
            os.makedirs(store_TFfam_dir)

        store_motif_common_name_dic = {}
        with open (ipt_motif_commonnm_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                motifID = col[0]
                motiffam = col[1]
                store_motif_common_name_dic[motifID] = motiffam

        store_target_TF_dic = {}
        store_target_motif_dic = {}
        with open (ipt_TF_motif_corr_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                motifID = col[0]
                TFID = col[1]

                if motifID in store_motif_common_name_dic:

                    motiffam = store_motif_common_name_dic[motifID]
                    if motiffam == eachfam:
                        store_target_motif_dic[motifID] = 1
                        store_target_TF_dic[TFID] = 1

        ##for the TF
        store_TF_dir = store_TFfam_dir + '/store_TF_dir'
        if not os.path.exists(store_TF_dir):
            os.makedirs(store_TF_dir)

        store_final_line_list = []
        first_line = 'name' + '\t' + 'geneID' + '\t' + 'type' + '\t' + 'tissue' + '\t' + 'common'
        store_final_line_list.append(first_line)

        for eachgene in store_target_TF_dic:
            final_line = eachgene + '\t' + eachgene + '\t' + 'TF' + '\t' + 'all' + '\t' + eachgene
            store_final_line_list.append(final_line)

        with open(opt_dir + '/opt_target_TF_markers.txt', 'w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        plot_TFs_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/plot_TFs_motifs/plot_genes_112723.R'

        cmd = 'Rscript ' + plot_TFs_script + \
              ' ' + ipt_gene_impute_acc_rds_fl + \
              ' ' + GAobj_rds_fl + \
              ' ' + opt_dir + '/opt_target_TF_markers.txt' + \
              ' ' + store_TF_dir + \
              ' ' + s0_s5_sub_lim
        print(cmd)
        subprocess.call(cmd, shell=True)

        ##for the motif
        store_motif_dir = store_TFfam_dir + '/store_motif_dir'
        if not os.path.exists(store_motif_dir):
            os.makedirs(store_motif_dir)

        store_final_line_list = []
        first_line = 'name' + '\t' + 'geneID' + '\t' + 'type' + '\t' + 'tissue' + '\t' + 'common'
        store_final_line_list.append(first_line)

        for eachmotif in store_target_motif_dic:
            final_line = eachmotif + '\t' + eachmotif + '\t' + 'TF' + '\t' + 'all' + '\t' + eachmotif
            store_final_line_list.append(final_line)

        with open(opt_dir + '/opt_target_motif_markers.txt', 'w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        plot_motifs_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/plot_TFs_motifs/plot_motifs_112723.R'

        cmd = 'Rscript ' + plot_motifs_script + \
              ' ' + ipt_motif_impute_acc_threesparse_fl + \
              ' ' + ipt_meta_fl + \
              ' ' + opt_dir + '/opt_target_motif_markers.txt' + \
              ' ' + store_motif_dir + \
              ' ' + s0_s5_sub_lim
        print(cmd)
        subprocess.call(cmd, shell=True)



##updating 120923
##select the potential TF and motif to be plotted based on correlation files and shared
##first checked the shared WRKY in the correlation files
##find the highest correlation to plot

def subfunction_select_target_TF_motif_in_all_spe_to_plot (ipt_correlation_dir,ipt_motif_commonnm_fl,ipt_target_plot_spe_str,ipt_target_plot_corr_cutoff,
                                                opt_dir):

    ##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_cross_species_091823/output_dir_103123/step00_basic_characters_dir/s0_s5_open_motifs_TFs_dir/store_corr_motif_TF_dir

    ipt_target_plot_spe_list = ipt_target_plot_spe_str.split(',')
    store_motif_spe_corr_dic = {}
    store_motif_gene_dic = {}
    for eachspe in ipt_target_plot_spe_list:

        ipt_spe_corr_fl = ipt_correlation_dir + '/' + eachspe + '/opt_corresponding_motif_TF.txt'

        with open (ipt_spe_corr_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                motifnm = col[0]
                genenm = col[1]
                corr = col[2]

                if float(corr) > float(ipt_target_plot_corr_cutoff):

                    if motifnm in store_motif_spe_corr_dic:
                        store_motif_spe_corr_dic[motifnm][eachspe] = 1
                    else:
                        store_motif_spe_corr_dic[motifnm] = {}
                        store_motif_spe_corr_dic[motifnm][eachspe] = 1

                    if motifnm in store_motif_gene_dic:
                        store_motif_gene_dic[motifnm][genenm] = 1
                    else:
                        store_motif_gene_dic[motifnm] = {}
                        store_motif_gene_dic[motifnm][genenm] = 1


    store_motif_common_name_dic = {}
    with open(ipt_motif_commonnm_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            motifID = col[0]
            motiffam = col[1]
            store_motif_common_name_dic[motifID] = motiffam

    store_final_line_list = []
    for eachmotifnm in store_motif_spe_corr_dic:

        if len(list(store_motif_spe_corr_dic[eachmotifnm].keys())) == len(ipt_target_plot_spe_list):

            genename_str = ','.join(list(store_motif_gene_dic[eachmotifnm].keys()))

            final_line = eachmotifnm + '\t' + genename_str + '\t' + store_motif_common_name_dic[eachmotifnm] + '\t' + eachmotifnm
            store_final_line_list.append(final_line)


    with open (opt_dir + '/opt_select_motif_TF_high_corr_' + ipt_target_plot_corr_cutoff + '_in_all_spe.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



##updating 123023
##we will plot the target plot only for the Uf
def subfunction_plot_motifs_for_Uf (input_required_scripts_dir,ipt_motif_TF_str,ipt_motif_impute_acc_threesparse_fl, ipt_meta_fl,
                        ipt_gene_impute_acc_rds_fl,GAobj_rds_fl,
                          opt_dir,s0_s5_sub_lim):


    ##we will check motif TF
    ##motif__TF,motif__TF
    ipt_motif_TF_list = ipt_motif_TF_str.split(',')


    ##opt_corresponding_motif_TF.txt is the ipt_TF_motif_corr_fl
    ##s0_s5_sub_plot_topnum_pairs is the number of
    #store_TF_motif_pair_coef_dic = {}
    #with open (ipt_TF_motif_corr_fl,'r') as ipt:
    #    for eachline in ipt:
    #        eachline = eachline.strip('\n')
    #        col = eachline.strip().split()
    #        motif_TF = col[0] + '__' + col[1]
    #        coef = float(col[2])
    #        store_TF_motif_pair_coef_dic[motif_TF] = coef

    ##select the top for building the candidate markers
    #res = dict(sorted(store_TF_motif_pair_coef_dic.items(), key=itemgetter(1), reverse=True)[:int(s0_s5_sub_plot_topnum_pairs)])

    store_target_motif_dic = {}
    for eachpair in ipt_motif_TF_list:
        mt = re.match('(.+)__(.+)',eachpair)
        motif = mt.group(1)
        gene = mt.group(2)
        store_target_motif_dic[motif] = 1

    store_final_line_list = []
    first_line = 'name' + '\t' + 'geneID' + '\t' +  'type' + '\t' + 'tissue' + '\t' + 'common'
    store_final_line_list.append(first_line)

    for eachmotif in store_target_motif_dic:
        final_line = eachmotif + '\t' + eachmotif + '\t' + 'TF' + '\t' + 'all' + '\t' + eachmotif
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_target_motif_markers.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    plot_motifs_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/plot_TFs_motifs/plot_motifs_112723.R'

    cmd = 'Rscript ' + plot_motifs_script + \
          ' ' + ipt_motif_impute_acc_threesparse_fl + \
          ' ' + ipt_meta_fl + \
          ' ' + opt_dir + '/opt_target_motif_markers.txt' + \
          ' ' + opt_dir + \
          ' ' + s0_s5_sub_lim
    print(cmd)
    subprocess.call(cmd,shell=True)




    store_target_gene_dic = {}
    for eachpair in ipt_motif_TF_list:
        mt = re.match('(.+)__(.+)', eachpair)
        motif = mt.group(1)
        gene = mt.group(2)
        store_target_gene_dic[gene] = 1

    store_final_line_list = []
    first_line = 'name' + '\t' + 'geneID' + '\t' + 'type' + '\t' + 'tissue' + '\t' + 'common'
    store_final_line_list.append(first_line)

    for eachgene in store_target_gene_dic:
        final_line = eachgene + '\t' + eachgene + '\t' + 'TF' + '\t' + 'all' + '\t' + eachgene
        store_final_line_list.append(final_line)

    with open(opt_dir + '/opt_target_TF_markers.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    plot_TFs_script = input_required_scripts_dir + '/s0_subfunctions_motifs_TFs/plot_TFs_motifs/plot_genes_112723.R'

    cmd = 'Rscript ' + plot_TFs_script + \
          ' ' + ipt_gene_impute_acc_rds_fl + \
          ' ' + GAobj_rds_fl + \
          ' ' + opt_dir + '/opt_target_TF_markers.txt' + \
          ' ' + opt_dir + \
          ' ' + s0_s5_sub_lim
    print(cmd)
    subprocess.call(cmd, shell=True)

























