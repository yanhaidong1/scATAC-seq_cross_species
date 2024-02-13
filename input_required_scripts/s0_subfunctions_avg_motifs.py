#!/usr/bin/env python

##updating 111823 we will use this script to check the motif coverage surrounding the ACR


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



def generate_exclude_region (input_te_gff_fl,input_acr_fl,opt_dir):

    store_exclude_region_line_list = []
    with open (input_te_gff_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split()
                type = col[2]
                if 'Gm' in col[0]:
                    if col[0].startswith('Gm'):

                        if type == 'repeat_region':
                            final_line = col[0] + '\t' + col[3] + '\t' + col[4]
                            store_exclude_region_line_list.append(final_line)

                else:
                    if type == 'repeat_region':
                        final_line = col[0] + '\t' + col[3] + '\t' + col[4]
                        store_exclude_region_line_list.append(final_line)

    with open (input_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            final_line = col[0] + '\t' + col[1] + '\t' + col[2]
            store_exclude_region_line_list.append(final_line)

    store_final_line_list = []
    for eachline in store_exclude_region_line_list:
        col = eachline.strip().split()
        if int(col[1]) < int(col[2]):
            store_final_line_list.append(eachline)

    with open (opt_dir + '/opt_exclude_region.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

def run_bedtools_shuffle (input_acr_fl,input_genome_size_fl,exclude_region_fl,opt_dir,simulated_time):

    store_final_line_list = []
    with open(input_acr_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '\t' + col[1] + '\t' + col[2]
            store_final_line_list.append(acrloc)
    with open (opt_dir + '/temp_3col_acr.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    store_all_candidate_region_line_list = []
    for i in range(int(simulated_time)):

        ##updating 121222 do not allow the overlapping
        cmd = 'bedtools shuffle -i ' + opt_dir + '/temp_3col_acr.txt' + ' -g ' + input_genome_size_fl + ' -excl ' + exclude_region_fl + ' -noOverlapping > ' \
              + opt_dir + '/temp_candidate_control_regions_' + str(i) + '.bed'
        #print(cmd)
        #subprocess.call(cmd,shell=True)

        cmd = 'bedtools shuffle -i ' + opt_dir + '/temp_3col_acr.txt' + ' -g ' + input_genome_size_fl + ' -excl ' + exclude_region_fl + ' > ' \
              + opt_dir + '/temp_candidate_control_regions_' + str(i) + '.bed'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##store all simulated control regions
        with open (opt_dir + '/temp_candidate_control_regions_' + str(i) + '.bed','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                eachline = eachline + '\t' + 'na' + '\t' + 'na'
                store_all_candidate_region_line_list.append(eachline)

    with open (opt_dir + '/temp_combine_candidate_regions.txt','+w') as opt:
        for eachline in store_all_candidate_region_line_list:
            opt.write(eachline + '\n')

    ##check whether the control bed file has duplicated regions
    store_acr_dup_dic = {}
    with open (opt_dir + '/temp_combine_candidate_regions.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            if acrnm in store_acr_dup_dic:
                store_acr_dup_dic[acrnm] += 1
            else:
                store_acr_dup_dic[acrnm] = 1

    with open (opt_dir + '/temp_checked_dup_acr.txt','w+') as opt:
        for eachline in store_acr_dup_dic:
            if store_acr_dup_dic[eachline] > 1:
                opt.write(eachline + '\t' + 'na' + '\t' + 'na' + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_combine_candidate_regions.txt > ' + opt_dir + '/opt_candidate_control_regions.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)





########################
##step01: create control
def subfunction_s1_create_control (ipt_spe_te_gff_fl,ipt_spe_acr_fl,ipt_spe_genome_size_fl,opt_dir):

    simulated_time = 1
    generate_exclude_region(ipt_spe_te_gff_fl, ipt_spe_acr_fl, opt_dir)
    run_bedtools_shuffle(ipt_spe_acr_fl, ipt_spe_genome_size_fl, opt_dir + '/opt_exclude_region.txt', opt_dir, simulated_time)

#########################
#########################



########################
##for check avg function
def reformat_acr_extend (input_ori_true_acr_fl,input_range,input_output_dir):

    store_final_line_list = []
    with open (input_ori_true_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            middle_point_acr = int((int(col[1]) + int(col[2]))/2)

            if int(middle_point_acr) - int(input_range) - 1 > 0:

                acr_st = middle_point_acr - int(input_range) - 1
                acr_ed = middle_point_acr + int(input_range) + 1

                other_line_list = []
                for i in range(3,len(col)):
                    other_line_list.append(col[i])

                final_line = col[0] + '\t' + str(acr_st) + '\t' + str(acr_ed) + '\t' + '\t'.join(other_line_list)
                store_final_line_list.append(final_line)

    with open (input_output_dir + '/opt_modi_acr.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


def intersect (input_ACR_fl,input_motif_fl,input_output_dir):

    ##Here the ACR should be five col

    cmd = 'bedtools intersect -wa -wb -a ' + input_ACR_fl + ' -b ' + input_motif_fl + ' > ' + input_output_dir + '/opt_ACR_motif.bed'
    print(cmd)
    subprocess.call(cmd,shell=True)

def modify_nm (opt_ACR_motif_fl,input_output_dir):

    store_total_motif_nm_dic = {}
    store_final_line_list = []
    with open (opt_ACR_motif_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            motif_nm = col[8]

            mt = re.match('(.+)_.+',motif_nm)
            concise_motif = mt.group(1)

            final_line = col[0]
            for i in range(1,8):
                final_line = final_line + '\t' + col[i]

            final_line = final_line + '\t' + concise_motif + '\t' + col[9] + '\t' + col[10] + '\t' + col[11]
            store_final_line_list.append(final_line)

            store_total_motif_nm_dic[concise_motif] = 1

    with open (input_output_dir + '/opt_ACR_motif_modinm.bed','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    with open (input_output_dir + '/opt_total_motif_nm.txt','w+') as opt:
        for eachline in store_total_motif_nm_dic:
            opt.write(eachline + '\n')




def reformat_acr_final_step (input_ori_true_acr_fl,input_output_dir):

    store_final_line_list = []
    with open (input_ori_true_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            arcnm = col[0] + '_' + col[1] + '_' + col[2]

            final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + arcnm + '\t' + col[5] + '\t' + col[6] + '\t' + col[7] +'\t' + col[8] + '\t' + col[9]  + '\t' + col[10]
            store_final_line_list.append(final_line)

    with open (input_output_dir + '/opt_modi_acr_motif.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



def subfunction_s2_check_avg_for_true_acr (opt_dir,ipt_spe_acr_fl,ipt_spe_motif_fl,ipt_convert2matrix_script,
                                        s0_s4_extend_acr_range_bp):

    ###########################
    ##true acr interaction file
    ##we first reformat the input_ACR_fl to check the up and down 2 kb information

    ##build the five col acr
    store_final_line_list = []
    with open (ipt_spe_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + 'na' + '\t' + 'na'
            store_final_line_list.append(acrloc)

    with open (opt_dir + '/temp_5col_acr.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    reformat_acr_extend(opt_dir + '/temp_5col_acr.txt', s0_s4_extend_acr_range_bp, opt_dir)
    intersect(opt_dir + '/opt_modi_acr.txt', ipt_spe_motif_fl, opt_dir)
    modify_nm(opt_dir + '/opt_ACR_motif.bed', opt_dir)
    reformat_acr_final_step(opt_dir + '/opt_ACR_motif_modinm.bed', opt_dir)

    ##for the true acr
    ipt_true_acr_motif_intersect_fl = opt_dir + '/opt_modi_acr_motif.txt'

    cmd = 'perl ' + ipt_convert2matrix_script + ' ' + ipt_true_acr_motif_intersect_fl + ' > ' + opt_dir + '/opt_position_count_mtx_true_acr.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)


def subfunction_s2_check_avg_for_control_acr (opt_dir,ipt_control_acr_fl,ipt_spe_motif_fl,s0_s4_extend_acr_range_bp,
                                           ipt_convert2matrix_script):

    reformat_acr_extend(ipt_control_acr_fl, s0_s4_extend_acr_range_bp, opt_dir)
    intersect(opt_dir + '/opt_modi_acr.txt', ipt_spe_motif_fl, opt_dir)
    modify_nm(opt_dir + '/opt_ACR_motif.bed', opt_dir)
    reformat_acr_final_step(opt_dir + '/opt_ACR_motif_modinm.bed', opt_dir)

    ##for the control acr
    ipt_control_acr_motif_intersect_fl = opt_dir + '/opt_modi_acr_motif.txt'

    cmd = 'perl ' + ipt_convert2matrix_script + ' ' + ipt_control_acr_motif_intersect_fl + ' > ' + opt_dir + '/opt_position_count_mtx_control_acr.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)


########################
##build a file for the R
def subfunction_s3_build_posi_coverage (opt_dir,ipt_position_count_mtx_true_acr_fl,
                                        ipt_position_count_mtx_control_acr_fl,
                                        ipt_R_script_for_build_coverage):

    #############
    ##s1 we first check number of true acr
    true_acr_count = 0
    with open(ipt_position_count_mtx_true_acr_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            true_acr_count += 1

    store_simu_acr_line_list = []
    control_acr_count = 0
    colnum = 0
    with open(ipt_position_count_mtx_control_acr_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            control_acr_count += 1
            col = eachline.strip().split()
            colnum = len(col)
            store_simu_acr_line_list.append(eachline)

    diff_acr_num = true_acr_count - control_acr_count

    for i in range(diff_acr_num):

        line_list = []
        for j in range(colnum):
            line_list.append('0')

        final_line = '\t'.join(line_list)
        store_simu_acr_line_list.append(final_line)

    with open(opt_dir + '/opt_position_count_mtx_control_acr_addequalACR.txt', 'w+') as opt:
        for eachline in store_simu_acr_line_list:
            opt.write(eachline + '\n')

    ########
    ##s2 use a R script to generate a dataframe to record the coverage
    ##s3_build_coverage_data_dir
    ##for the true
    ipt_script = ipt_R_script_for_build_coverage

    cmd = 'Rscript ' + ipt_script + \
          ' ' + ipt_position_count_mtx_true_acr_fl + \
          ' ' + opt_dir + \
          ' ' + 'true'
    subprocess.call(cmd, shell=True)

    cmd = 'Rscript ' + ipt_script + \
          ' ' + opt_dir + '/opt_position_count_mtx_control_acr_addequalACR.txt' + \
          ' ' + opt_dir + \
          ' ' + 'control'
    subprocess.call(cmd, shell=True)


















