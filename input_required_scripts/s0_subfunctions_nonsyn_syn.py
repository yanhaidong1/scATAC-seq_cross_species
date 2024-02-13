#!/usr/bin/env python

##this script is to check the proportion of syntenic and nonsyntenic ACRs of the cell-type-specific ones
##and do the enrichment test


import re
import glob
import sys
import subprocess
import os
from statistics import mean
import scipy.stats as stats
import random


##input_rice_acr_add_celltype_fl = sys.argv[10]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_6_characterize_ACRs_060322/output_dir_v8.2_update_051323/step01_ACR_summary/s1_loc_analysis/opt_ACR_sorted_rmdup_addallcelltype.txt



def create_bins(lower_bound, upper_bound, width):
    """ create_bins returns an equal-width (distance) partitioning.
        It returns an ascending list of tuples, representing the intervals.
        A tuple bins[i], i.e. (bins[i][0], bins[i][1])  with i > 0
        and i < quantity, satisfies the following conditions:
            (1) bins[i][0] + width == bins[i][1]
            (2) bins[i-1][0] + width == bins[i][0] and
                bins[i-1][1] + width == bins[i][1]
    """
    bins = []
    current = lower_bound
    while current <= upper_bound:
        bins.append((current, current + width + 1))
        current += width

    return bins

def subfunctions_check_syn_nonsyn_celltypebin_prop (opt_dir,input_rice_acr_add_celltype_fl,ipt_all_ACR_dic,ipt_target_ACR_dic,s0_s1_width_bin,
                                                target_prefix):

    ##target ACR dic could be the ACR within the syntenic or ACRs not within the syntenic
    all_ACR_list = list(ipt_all_ACR_dic.keys())

    ACR_contain_CNS_num = len(list(ipt_target_ACR_dic.keys()))

    ##Here we will check proportion of cell types
    cell_num_list = []
    store_ACR_celltypestr_dic = {}
    store_ACR_celltypenum_dic = {}
    store_allcelltype_dic = {}
    with open(input_rice_acr_add_celltype_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            celltypenum = col[5]
            if celltypenum != 'none':
                cell_num_list.append(int(celltypenum))
                store_ACR_celltypestr_dic[acrnm] = col[6]
                store_ACR_celltypenum_dic[acrnm] = int(celltypenum)

                all_celltype_list = col[6].split(',')
                for eachcelltype in all_celltype_list:
                    store_allcelltype_dic[eachcelltype] = 1

    max_cell_num = max(cell_num_list)

    lower_bound = 0
    upper_bound = max_cell_num
    width = int(s0_s1_width_bin)

    bins_final = create_bins(lower_bound, upper_bound, width)
    print(bins_final)

    store_bin_simulate_prop_list_dic = {}
    store_bin_simulate_eachcelltypePropDic_list_dic = {}
    ##Here we just randomly check 5 time to get the real other than the control
    for i in range(int(2)):

        store_rd_ACR_list = random.sample(all_ACR_list, ACR_contain_CNS_num)

        ##for each bin
        for eachbin in bins_final:

            binID = '[' + str(eachbin[0] + 1) + ',' + str(eachbin[1] - 1) + ']'

            bin_st = eachbin[0]
            bin_ed = eachbin[1]

            ##extract the ACR based on the bine
            store_rd_target_bin_ACR_list = []

            for eachACR in store_rd_ACR_list:
                if eachACR in store_ACR_celltypenum_dic:
                    celltypenum = store_ACR_celltypenum_dic[eachACR]
                    if celltypenum > bin_st and celltypenum < bin_ed:
                        store_rd_target_bin_ACR_list.append(eachACR)

            store_rd_bin_ACR_num = len(store_rd_target_bin_ACR_list)

            prop = store_rd_bin_ACR_num / ACR_contain_CNS_num

            if binID in store_bin_simulate_prop_list_dic:
                store_bin_simulate_prop_list_dic[binID].append(prop)
            else:
                store_bin_simulate_prop_list_dic[binID] = []
                store_bin_simulate_prop_list_dic[binID].append(prop)

            #################
            ##updating 101623 we will generate a prop for each number of cell type
            store_simulate_ACR_celltypecomp_dic = {}
            for eachACR in store_rd_ACR_list:
                if eachACR in store_ACR_celltypenum_dic:
                    celltypenum = store_ACR_celltypenum_dic[eachACR]
                    if celltypenum > bin_st and celltypenum < bin_ed:

                        celltypestr = store_ACR_celltypestr_dic[eachACR]
                        celltypelist = celltypestr.split(',')
                        for eachcelltype in celltypelist:
                            if eachcelltype in store_simulate_ACR_celltypecomp_dic:
                                store_simulate_ACR_celltypecomp_dic[eachcelltype] += 1
                            else:
                                store_simulate_ACR_celltypecomp_dic[eachcelltype] = 1

            ##Here we will use the rd bin ACR to be the total
            store_celltype_prop_dic = {}
            for eachcelltype in store_simulate_ACR_celltypecomp_dic:
                celltypenum = store_simulate_ACR_celltypecomp_dic[eachcelltype]

                ##Here we will use the contain CNS num
                prop = celltypenum / ACR_contain_CNS_num

                # prop = celltypenum/store_rd_bin_ACR_num
                store_celltype_prop_dic[eachcelltype] = prop

            if binID in store_bin_simulate_eachcelltypePropDic_list_dic:
                store_bin_simulate_eachcelltypePropDic_list_dic[binID].append(store_celltype_prop_dic)
            else:
                store_bin_simulate_eachcelltypePropDic_list_dic[binID] = []
                store_bin_simulate_eachcelltypePropDic_list_dic[binID].append(store_celltype_prop_dic)

    ##this file is right
    print(store_bin_simulate_eachcelltypePropDic_list_dic)

    ##check the real number
    store_bin_ACRinCNS_num_dic = {}
    store_bin_real_celltypeNumDic_dic = {}
    for eachbin in bins_final:
        binID = '[' + str(eachbin[0] + 1) + ',' + str(eachbin[1] - 1) + ']'

        bin_st = eachbin[0]
        bin_ed = eachbin[1]

        store_CNS_bin_ACR_list = []
        for eachACR in ipt_target_ACR_dic:
            if eachACR in store_ACR_celltypenum_dic:
                celltypenum = store_ACR_celltypenum_dic[eachACR]
                if celltypenum > bin_st and celltypenum < bin_ed:
                    store_CNS_bin_ACR_list.append(eachACR)

        ACR_in_CNS_num = len(store_CNS_bin_ACR_list)
        store_bin_ACRinCNS_num_dic[binID] = ACR_in_CNS_num

        #################
        ##updating 101623 we will generate a prop for each cell type
        store_real_ACR_celltypecomp_dic = {}
        for eachACR in store_CNS_bin_ACR_list:

            celltypestr = store_ACR_celltypestr_dic[eachACR]
            celltypelist = celltypestr.split(',')
            for eachcelltype in celltypelist:
                if eachcelltype in store_real_ACR_celltypecomp_dic:
                    store_real_ACR_celltypecomp_dic[eachcelltype] += 1
                else:
                    store_real_ACR_celltypecomp_dic[eachcelltype] = 1

        ##Here we will use the rd bin ACR to be the total
        store_celltype_num_dic = {}
        for eachcelltype in store_real_ACR_celltypecomp_dic:
            celltypenum = store_real_ACR_celltypecomp_dic[eachcelltype]
            store_celltype_num_dic[eachcelltype] = celltypenum

        store_bin_real_celltypeNumDic_dic[binID] = store_celltype_num_dic

    ##check for eachbin
    store_final_line_Rplot_list = []
    for eachbinID in store_bin_simulate_prop_list_dic:
        prop_list = store_bin_simulate_prop_list_dic[eachbinID]

        for eachprop in prop_list:
            final_line_plotR = eachbinID + '\t' + 'control' + '\t' + str(eachprop)
            store_final_line_Rplot_list.append(final_line_plotR)

        real_num = store_bin_ACRinCNS_num_dic[eachbinID]

        final_line_plotR = eachbinID + '\t' + 'real' + '\t' + str(real_num / ACR_contain_CNS_num)
        store_final_line_Rplot_list.append(final_line_plotR)

    with open(opt_dir + '/opt_' + target_prefix + '_in_diffbins_plot_R.txt',
              'w+') as opt:
        for eachline in store_final_line_Rplot_list:
            opt.write(eachline + '\n')


##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_cross_species_091823/output_dir_103123/step00_basic_characters_dir/s0_s1_open_check_syntenic_block_dir/opt_allrice_syntenic_region_add_distance_sorted.txt

def subfunction_check_overlap_number_ACRs_in_syn_non_syn(input_rice_acr_add_celltype_fl,input_rice_atlas_acr_fl,ipt_syntenic_block_fl, opt_dir,s0_s1_width_bin):
    store_final_line_dic = {}
    with open(ipt_syntenic_block_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            final_line = col[0] + '\t' + col[1] + '\t' + col[2]
            store_final_line_dic[final_line] = 1

    with open(opt_dir + '/temp_rice_syntenic_nosort.txt', 'w+') as opt:
        for eachline in store_final_line_dic:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_rice_syntenic_nosort.txt > ' + opt_dir + '/temp_rice_syntenic_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = 'cut -f 1-3 ' + input_rice_acr_add_celltype_fl + ' > ' + opt_dir + '/temp_rice_atlas_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_rice_atlas_acr.txt > ' + opt_dir + '/temp_rice_atlas_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_rice_atlas_acr_sorted.txt' + ' -b ' + opt_dir + '/temp_rice_syntenic_sorted.txt > ' + \
          opt_dir + '/temp_intersect_rice_ACR_syntenic.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##store the ACR within CNS
    store_ACR_within_syntenic_dic = {}
    with open(opt_dir + '/temp_intersect_rice_ACR_syntenic.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            store_ACR_within_syntenic_dic[acrnm] = 1

    store_all_ACR_dic = {}
    store_ACR_not_within_syntenic_dic = {}
    with open (opt_dir + '/temp_rice_atlas_acr_sorted.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            store_all_ACR_dic[acrloc] = 1
            if acrloc not in store_ACR_within_syntenic_dic:
                store_ACR_not_within_syntenic_dic[acrloc] = 1

    ##we will check the syntenic ones
    target_prefix = 'Syn'
    subfunctions_check_syn_nonsyn_celltypebin_prop(opt_dir, input_rice_acr_add_celltype_fl, store_all_ACR_dic,
                                                   store_ACR_within_syntenic_dic, s0_s1_width_bin,
                                                   target_prefix)

    target_prefix = 'NoSyn'
    subfunctions_check_syn_nonsyn_celltypebin_prop(opt_dir, input_rice_acr_add_celltype_fl, store_all_ACR_dic,
                                                   store_ACR_not_within_syntenic_dic, s0_s1_width_bin,
                                                   target_prefix)

    store_final_line_list = []
    with open (input_rice_atlas_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            if acrloc in store_ACR_within_syntenic_dic:
                final_cate = 'Syn'
            else:
                final_cate = 'NonSyn'
            final_line = eachline + '\t' + final_cate
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_rice_atlas_acr_add_syn_or_not.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


##updating 021224
##we will check the size of syntenic regions
def subfunction_check_size_of_syntenic_region (opt_dir):

    ipt_region_fl = opt_dir + '/temp_rice_syntenic_sorted.txt'

    cmd = 'bedtools merge -i ' + ipt_region_fl + ' > ' + opt_dir + '/temp_rice_syntenic_sorted_merge.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    total_len = 0
    with open (opt_dir + '/temp_rice_syntenic_sorted_merge.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            diff = int(col[2]) - int(col[1])
            total_len += diff


    with open (opt_dir + '/temp_rice_syntenic_region_total_len.txt','w+') as opt:
        opt.write('total len ' + str(total_len) + '\n')



























