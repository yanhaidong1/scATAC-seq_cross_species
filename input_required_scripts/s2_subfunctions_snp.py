#!/usr/bin/env python


##updating 122723 this is for the snp other than the cns
##updating 120423 we will do the enrichment test for the restricted ACRs
##updating 112523 we will check the broad and CT peaks overlapping with 128 cell types
##updating 111723 we will add a new function to check the accessiblility of the rice altas peaks for the syntenic peaks
##updating 111523 we will add the categories of ACR nearby genes
##updating 111523 check if broad ACRs were relative higher enriched compared to the others


import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats
import random
from statistics import mean

##updating 122723
def subfunction_add_snp_num_to_acr (ipt_snp_fl,ipt_final_summary_fl,spe1_prefix,spe2_prefix,opt_dir):

    ##build the snp bed fl
    store_final_line_list = []
    with open (ipt_snp_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if not eachline.startswith('#'):

                chrnm = 'Chr' + col[0]
                pos = col[1]

                final_line = chrnm + '\t' + pos + '\t' + str(int(pos) + 1)
                store_final_line_list.append(final_line)

    with open (opt_dir + '/temp_snp_posi_fl.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_snp_posi_fl.txt > ' + opt_dir + '/temp_snp_posi_fl_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)


    store_acr_line_list = []
    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                acr_line = '\t'.join(col[1].split('_'))
                store_acr_line_list.append(acr_line)

    with open (opt_dir + '/temp_acr.txt','w+') as opt:
        for eachline in store_acr_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_acr.txt > ' + opt_dir + '/temp_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_acr_sorted.txt' + \
          ' -b ' + opt_dir + '/temp_snp_posi_fl_sorted.txt > ' + opt_dir + '/temp_intersect_acr_snp.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_acr_snp_num_dic = {}
    with open (opt_dir + '/temp_intersect_acr_snp.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acr = col[0] + '_' + col[1] + '_' + col[2]
            snp = col[3] + '_' + col[4] + '_' + col[5]

            if acr in store_acr_snp_num_dic:
                store_acr_snp_num_dic[acr][snp] = 1
            else:
                store_acr_snp_num_dic[acr] = {}
                store_acr_snp_num_dic[acr][snp] = 1

    store_final_line_list = []
    count = 0

    with open(ipt_final_summary_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                ##store the different cate of acr
                spe1_acrID = col[0]
                spe1_acrloc = col[1]

                if spe1_acrloc in store_acr_snp_num_dic:
                    snp_num = len(list(store_acr_snp_num_dic[spe1_acrloc].keys()))
                else:
                    snp_num = 0

                final_line = eachline + '\t' + str(snp_num)

            else:
                final_line = eachline + '\t' + spe1_prefix + '_captureSNPnum'

            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_to_' + spe2_prefix + '_summary_add_' + spe1_prefix + '_snp_num.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




