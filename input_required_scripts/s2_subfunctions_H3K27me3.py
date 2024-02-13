#!/usr/bin/env python

##updating 122723 we wil
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




def subfunction_plot_H3K27me3_on_species_specific_ACRs (ipt_final_summary_spe_fl_dic,ipt_H3K27me3_bw_fl,opt_dir,input_core_num):


    store_final_acr_fl_list = []
    for eachspe2 in ipt_final_summary_spe_fl_dic:

        ipt_final_summary_fl = ipt_final_summary_spe_fl_dic[eachspe2]
        ##build the specific ACRs
        store_final_line_list = []
        count = 0
        with open (ipt_final_summary_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1
                if count != 1:
                    spe2_region = col[4]
                    if spe2_region == 'none':
                        spe1_acrloc = col[1]
                        spe1_loc = '\t'.join(spe1_acrloc.split('_'))
                        store_final_line_list.append(spe1_loc)

        with open (opt_dir + '/temp_acr.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_acr.txt > ' + opt_dir + '/' + eachspe2
        print(cmd)
        subprocess.call(cmd,shell=True)

        store_final_acr_fl_list.append(opt_dir + '/' + eachspe2)


    cmd = 'computeMatrix reference-point' + \
          ' --referencePoint center' + \
          ' -a 2000 -b 2000' + \
          ' -S ' + ipt_H3K27me3_bw_fl + \
          ' -R ' + ' '.join(store_final_acr_fl_list) + \
          ' -p ' + input_core_num + \
          ' --missingDataAsZero' + \
          ' -out ' + opt_dir + '/temp_H3K27m3_allspe_specific.tab.gz'
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = 'plotHeatmap' + \
          ' -m ' + opt_dir + '/temp_H3K27m3_allspe_specific.tab.gz' + \
          ' --perGroup' + \
          ' -out ' + opt_dir + '/opt_H3K27m3_allspe_species_specific.png'
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = 'plotHeatmap' + \
          ' -m ' + opt_dir + '/temp_H3K27m3_allspe_specific.tab.gz' + \
          ' --perGroup' + \
          ' --plotFileFormat pdf' + \
          ' -out ' + opt_dir + '/opt_H3K27m3_allspe_species_specific.pdf'
    print(cmd)
    subprocess.call(cmd, shell=True)







