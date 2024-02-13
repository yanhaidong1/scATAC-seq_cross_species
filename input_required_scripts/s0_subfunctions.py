#!/usr/bin/env python

##updating 111623 we will add a function to check the TE
##updating 111623 we will add a total of TE

import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats

def subfunction_check_syntenic_blocks (ipt_spe1_all_os_ACRs_fl, ipt_spe2_os_ACRs_fl,input_spe1_H3K27me3_fl,opt_dir,spe1_prefix,spe2_prefix):

    ##We aim to use the ipt_spe2_os_ACRs_fl to find the region
    store_spe2_region_dic = {}
    with open (ipt_spe2_os_ACRs_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            spe2_syn_region = col[3]
            store_spe2_region_dic[spe2_syn_region] = 1

    store_final_line_list = []
    store_region_ID_gene_pair_sted_dic = {}
    with open (ipt_spe1_all_os_ACRs_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            syntenic_region_ID = col[-1]
            if syntenic_region_ID in store_spe2_region_dic:
                store_final_line_list.append(eachline)

                gene_line = col[0] + '_' + col[1] + '_' + col[2] + '_' + col[3]

                if syntenic_region_ID in store_region_ID_gene_pair_sted_dic:
                    store_region_ID_gene_pair_sted_dic[syntenic_region_ID][gene_line] = 1
                    #store_region_ID_gene_pair_sted_dic[syntenic_region_ID].append(gene_line)
                else:

                    store_region_ID_gene_pair_sted_dic[syntenic_region_ID] = {}
                    store_region_ID_gene_pair_sted_dic[syntenic_region_ID][gene_line] = 1
                    #store_region_ID_gene_pair_sted_dic[syntenic_region_ID].append(gene_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_genes.all_os_ACRs.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##check the distance
    store_bug_line_list = []
    store_final_line_list = []
    for eachsyntenic_region_ID in store_region_ID_gene_pair_sted_dic:

        gene_line_dic = store_region_ID_gene_pair_sted_dic[eachsyntenic_region_ID]
        if len(list(gene_line_dic.keys())) == 2:

            store_location_list = []
            gene_chr = ''
            for eachgeneline in list(gene_line_dic.keys()):
                gene_col = eachgeneline.split('_')
                gene_st = gene_col[1]
                gene_ed = gene_col[2]
                store_location_list.append(int(gene_st))
                store_location_list.append(int(gene_ed))
                gene_chr = gene_col[0]

            store_location_list.sort()

            loc_2nd = store_location_list[1]
            loc_3rd = store_location_list[2]
            loc_1st = store_location_list[0]
            loc_4th = store_location_list[3]

            ##updating 102623 we will use the 1st and 4th
            #distance = abs(loc_3rd - loc_2nd)
            distance = abs(loc_4th - loc_1st)
            final_line = gene_chr + '\t' + str(loc_1st) + '\t' + str(loc_4th) + '\t' + eachsyntenic_region_ID + '\t' + str(distance)
            store_final_line_list.append(final_line)

        else:
            store_bug_line_list.append(eachsyntenic_region_ID + '\t' + '\t'.join(list(gene_line_dic.keys())))

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_bug_line.txt', 'w+') as opt:
        for eachline in store_bug_line_list:
            opt.write(eachline + '\n')

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_add_distance.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_add_distance.txt > ' + \
          opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_add_distance_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##check the number of ACRs within the gene pair
    cmd = 'cut -f 1-3 ' + input_spe1_H3K27me3_fl + ' > ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt > ' + opt_dir + '/temp_' + spe1_prefix + '_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + spe1_prefix + '_acr_sorted.txt' + \
          ' -b ' + opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_add_distance_sorted.txt' + \
          ' > ' + opt_dir + '/temp_intersect_acr_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_syntenic_region_spe1_ACR_dic = {}
    store_intersected_acr_dic = {}
    with open (opt_dir + '/temp_intersect_acr_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            ACRnm = col[0] + '_' + col[1] + '_' + col[2]
            syntenic_regionID = col[6]
            store_intersected_acr_dic[ACRnm] = 1

            if syntenic_regionID in store_syntenic_region_spe1_ACR_dic:
                store_syntenic_region_spe1_ACR_dic[syntenic_regionID][ACRnm] = 1
            else:
                store_syntenic_region_spe1_ACR_dic[syntenic_regionID] = {}
                store_syntenic_region_spe1_ACR_dic[syntenic_regionID][ACRnm] = 1

    with open (opt_dir + '/opt_ACRs_in_syntenic_regions.txt','w+') as opt:
        for eachacr in store_intersected_acr_dic:
            opt.write(eachacr + '\n')

    store_final_line_list = []
    for eachregion in store_syntenic_region_spe1_ACR_dic:
        ACR_str = ','.join(list(store_syntenic_region_spe1_ACR_dic[eachregion].keys()))
        ACR_num = len(list(store_syntenic_region_spe1_ACR_dic[eachregion].keys()))
        final_line = eachregion + '\t' + ACR_str + '\t' + str(ACR_num)
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_contain_ACRs.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

##updating 102023
def subfunction_build_syntenic_bed (ipt_spe1_all_os_ACRs_fl,input_spe1_H3K27me3_fl,opt_dir,spe1_prefix):

    store_final_line_list = []
    store_region_ID_gene_pair_sted_dic = {}
    with open(ipt_spe1_all_os_ACRs_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            syntenic_region_ID = col[-1]

            store_final_line_list.append(eachline)

            gene_line = col[0] + '_' + col[1] + '_' + col[2] + '_' + col[3]

            if syntenic_region_ID in store_region_ID_gene_pair_sted_dic:
                store_region_ID_gene_pair_sted_dic[syntenic_region_ID][gene_line] = 1
                # store_region_ID_gene_pair_sted_dic[syntenic_region_ID].append(gene_line)
            else:
                store_region_ID_gene_pair_sted_dic[syntenic_region_ID] = {}
                store_region_ID_gene_pair_sted_dic[syntenic_region_ID][gene_line] = 1
                # store_region_ID_gene_pair_sted_dic[syntenic_region_ID].append(gene_line)


    ##check the distance
    store_final_line_list = []
    for eachsyntenic_region_ID in store_region_ID_gene_pair_sted_dic:

        gene_line_dic = store_region_ID_gene_pair_sted_dic[eachsyntenic_region_ID]
        if len(list(gene_line_dic.keys())) == 2:

            store_location_list = []
            gene_chr = ''
            for eachgeneline in list(gene_line_dic.keys()):
                gene_col = eachgeneline.split('_')
                gene_st = gene_col[1]
                gene_ed = gene_col[2]
                store_location_list.append(int(gene_st))
                store_location_list.append(int(gene_ed))
                gene_chr = gene_col[0]

            store_location_list.sort()

            loc_2nd = store_location_list[1]
            loc_3rd = store_location_list[2]
            loc_1st = store_location_list[0]
            loc_4th = store_location_list[3]

            distance = abs(loc_3rd - loc_2nd)
            final_line = gene_chr + '\t' + str(loc_1st) + '\t' + str(
                loc_4th) + '\t' + eachsyntenic_region_ID + '\t' + str(distance)
            store_final_line_list.append(final_line)

        ##We do not consider this case, as there are three gene cases showing the genes overlapping with each other
        else:
            ##if we meet three or more
            ##we still need to find the location

            store_location_list = []
            gene_chr = ''
            for eachgeneline in list(gene_line_dic.keys()):
                gene_col = eachgeneline.split('_')
                gene_st = gene_col[1]
                gene_ed = gene_col[2]
                gene_chr = gene_col[0]

                store_location_list.append(int(gene_st))
                store_location_list.append(int(gene_ed))

            store_location_list.sort()

            total_site_len = len(store_location_list)
            gene_num = total_site_len/2

            print(store_location_list)
            print(gene_num)

            for i in range(int(gene_num) - 1):

                inter_st_posi_num = i*2 + 2
                inter_ed_posi_num = i*2 + 3
                print(inter_st_posi_num)
                print(inter_ed_posi_num)

                inter_loc_st = store_location_list[inter_st_posi_num - 1]
                inter_loc_ed = store_location_list[inter_ed_posi_num - 1]

                distance = abs(inter_loc_ed - inter_loc_st)

                ##
                #final_line = gene_chr + '\t' + str(inter_loc_st) + '\t' + str(inter_loc_ed) + '\t' + eachsyntenic_region_ID + '_' + str(i) + '\t' + str(distance)
                #store_final_line_list.append(final_line)

            ##1 2 3 4 5 6
            ##23 and 45

            ##12 34 56 78
            ##2 3 and 4 5 and 6 7


    with open(opt_dir + '/opt_all' + spe1_prefix +  '_syntenic_region_add_distance.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt_all' + spe1_prefix +  '_syntenic_region_add_distance.txt > ' + \
          opt_dir + '/opt_all' + spe1_prefix + '_syntenic_region_add_distance_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##intersect with the syntenic regions
    cmd = 'cut -f 1-3 ' + input_spe1_H3K27me3_fl + ' > ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt' + ' -b ' + \
          opt_dir + '/opt_all' + spe1_prefix + '_syntenic_region_add_distance_sorted.txt > ' + \
          opt_dir + '/temp_intersect_all_' + spe1_prefix + '_ACR_syntenic.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_peak_dic = {}
    with open (opt_dir + '/temp_intersect_all_' + spe1_prefix + '_ACR_syntenic.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            peaknm = col[0] + '_' + col[1] + '_' + col[2]
            store_peak_dic[peaknm] = 1

    with open (opt_dir + '/opt_' + spe1_prefix + '_peak_in_syntenic_region.txt','w+') as opt:
        for eachpeak in store_peak_dic:
            opt.write(eachpeak + '\n')


##updating 111623
def subfunction_summarize_te (store_type_dic,store_te_loc_nm_dic,store_te_prop_line_list):

    ##now we will check the proportion of TE
    #store_te_prop_line_list = []
    for eachdloc_type in store_type_dic:
        TE_loc_dic = store_type_dic[eachdloc_type]
        total_TE_num = len(list(TE_loc_dic.keys()))

        store_TEfam_dic = {}
        for eachteloc in TE_loc_dic:
            te_nm = store_te_loc_nm_dic[eachteloc]

            if te_nm in store_TEfam_dic:
                store_TEfam_dic[te_nm] += 1
            else:
                store_TEfam_dic[te_nm] = 1

        for eachte_allfam in store_TEfam_dic:
            te_num = store_TEfam_dic[eachte_allfam]
            te_prop = str(te_num / total_TE_num)

            if eachte_allfam != 'target_site_duplication' and \
                    eachte_allfam != 'long_terminal_repeat':

                new_te_all_fam_nm = ''
                if eachte_allfam == 'repeat_region':
                    new_te_all_fam_nm = 'Others'
                if eachte_allfam == 'helitron':
                    new_te_all_fam_nm = 'Helitron'
                if eachte_allfam == 'Mutator_TIR_transposon':
                    new_te_all_fam_nm = 'TIR_Mutator'
                if eachte_allfam == 'Tc1_Mariner_TIR_transposon':
                    new_te_all_fam_nm = 'TIR_Tc1_Mariner'
                if eachte_allfam == 'CACTA_TIR_transposon':
                    new_te_all_fam_nm = 'TIR_CACTA'
                if eachte_allfam == 'hAT_TIR_transposon':
                    new_te_all_fam_nm = 'TIR_hAT'
                if eachte_allfam == 'Gypsy_LTR_retrotransposon':
                    new_te_all_fam_nm = 'LTR_Cypsy'
                if eachte_allfam == 'PIF_Harbinger_TIR_transposon':
                    new_te_all_fam_nm = 'TIR_PIF_Harbinger'
                if eachte_allfam == 'LTR_retrotransposon':
                    new_te_all_fam_nm = 'LTR_Other'
                if eachte_allfam == 'Copia_LTR_retrotransposon':
                    new_te_all_fam_nm = 'LTR_Copia'
                if eachte_allfam == 'LINE_element':
                    new_te_all_fam_nm = 'LINE'
                if eachte_allfam == 'TRIM':
                    new_te_all_fam_nm = 'TRIM'
                if eachte_allfam == 'terminal_inverted_repeat_element':
                    new_te_all_fam_nm = 'TIRE'
                if eachte_allfam == 'centromeric_repeat':
                    new_te_all_fam_nm = 'CENTR'
                if eachte_allfam == 'non_LTR_retrotransposon':
                    new_te_all_fam_nm = 'nonLTRretroTE'
                if eachte_allfam == 'SINE_element':
                    new_te_all_fam_nm = 'SINE'

                if new_te_all_fam_nm == '':
                    new_te_all_fam_nm = eachte_allfam

                ##updating 113022
                ##add the super family infomration
                if 'DNA' in new_te_all_fam_nm:

                    if 'satellite_DNA_Satellite.rice' not in new_te_all_fam_nm:

                        if 'terminal_inverted_repeat_element_' not in new_te_all_fam_nm:

                            if 'helitron' in new_te_all_fam_nm:
                                new_super_fam_nm = 'Helitron'
                            else:
                                new_super_fam_nm = 'nMITE'

                        else:
                            new_super_fam_nm = 'nMITEothers_' + new_te_all_fam_nm

                    else:
                        new_super_fam_nm = 'SateliteDNA'

                else:

                    if 'MITE' in new_te_all_fam_nm:
                        new_super_fam_nm = 'MITE'

                    else:
                        if 'LTR' in new_te_all_fam_nm:

                            if 'non_LTR' not in new_te_all_fam_nm:
                                new_super_fam_nm = 'LTR'
                            else:
                                new_super_fam_nm = 'nLTRothers_' + new_te_all_fam_nm

                        else:

                            if 'LINE' in new_te_all_fam_nm:
                                new_super_fam_nm = 'LINE'
                            else:
                                if 'SINE' in new_te_all_fam_nm:
                                    new_super_fam_nm = 'SINE'

                                else:
                                    new_super_fam_nm = 'Others_' + new_te_all_fam_nm

                final_line =  eachdloc_type + '\t' + new_te_all_fam_nm + '\t' + str(
                    te_num) + '\t' + te_prop + '\t' + new_super_fam_nm
                store_te_prop_line_list.append(final_line)

def subfunction_check_TE_composition (ipt_spe_te_fl,ipt_acr_fl,spe_prefix,opt_dir):


    store_final_line = []
    with open(ipt_spe_te_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split()

                ##updating 070122
                annot_col = col[8].split(';')
                annot_dic = {}
                for eachannotstr in annot_col:
                    if '=' in eachannotstr:
                        mt = re.match('(.+)=(.+)', eachannotstr)
                        left = mt.group(1)
                        right = mt.group(2)
                        annot_dic[left] = right

                if 'Classification' in annot_dic:
                    subclass = annot_dic['Classification']
                else:
                    subclass = 'NoSub'

                new_subclass = subclass.replace('/', '.')

                final_line = col[0] + '\t' + col[3] + '\t' + col[4] + '\t' + col[6] + '\t' + col[2] + '_' + new_subclass
                store_final_line.append(final_line)

    with open(opt_dir + '/temp_TE.bed', 'w+') as opt:
        for eachline in store_final_line:
            opt.write(eachline + '\n')

    ##order these two file
    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_TE.bed > ' + opt_dir + '/temp_TE_sorted.bed'
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = 'sort -k1,1V -k2,2n ' + ipt_acr_fl + ' > ' + opt_dir + '/temp_ACR_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    # intersect with ACR
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_TE_sorted.bed -b ' + opt_dir + '/temp_ACR_sorted.txt' + ' > ' + \
          opt_dir + '/temp_TE_intersect_ACR.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    store_te_loc_nm_dic = {}
    with open(opt_dir + '/temp_TE_intersect_ACR.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            TE_loc = col[0] + '_' + col[1] + '_' + col[2]
            TE_nm = col[4]
            store_te_loc_nm_dic[TE_loc] = TE_nm

    store_type_dic = {}
    store_loctype_loc = {}
    with open(opt_dir + '/temp_TE_intersect_ACR.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')

            ##updating 111623
            ##check the loc_type
            loc_type = col[8]

            mt = re.match('.+;(.+)',loc_type)
            celltype = mt.group(1)
            if celltype == 'broadly_accessible':
                final_type = 'Broad'
            else:
                final_type = 'CT'

            TE_loc = col[0] + '_' + col[1] + '_' + col[2]

            if final_type in store_type_dic:
                store_type_dic[final_type][TE_loc] = 1
            else:
                store_type_dic[final_type] = {}
                store_type_dic[final_type][TE_loc] = 1

            te_acr_loc = col[5] + '_' + col[6] + '_' + col[7]
            store_loctype_loc[te_acr_loc] = final_type

    ##Here we will store the CT and all at same time
    store_te_prop_line_list = []
    subfunction_summarize_te(store_type_dic, store_te_loc_nm_dic, store_te_prop_line_list)

    ##store the whole genome te information
    store_final_line = []
    with open(ipt_spe_te_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split()

                ##updating 070122
                annot_col = col[8].split(';')
                annot_dic = {}
                for eachannotstr in annot_col:
                    if '=' in eachannotstr:
                        mt = re.match('(.+)=(.+)', eachannotstr)
                        left = mt.group(1)
                        right = mt.group(2)
                        annot_dic[left] = right

                if 'Classification' in annot_dic:
                    subclass = annot_dic['Classification']
                else:
                    subclass = 'NoSub'

                new_subclass = subclass.replace('/', '.')

                final_line = col[0] + '\t' + col[3] + '\t' + col[4] + '\t' + col[6] + '\t' + col[2] + '_' + new_subclass
                store_final_line.append(final_line)

    with open(opt_dir + '/temp_TE.bed', 'w+') as opt:
        for eachline in store_final_line:
            opt.write(eachline + '\n')

    store_te_loc_nm_dic = {}
    with open(opt_dir + '/temp_TE.bed', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            TE_loc = col[0] + '_' + col[1] + '_' + col[2]
            TE_nm = col[4]
            store_te_loc_nm_dic[TE_loc] = TE_nm

    store_type_dic = {}
    with open(opt_dir + '/temp_TE.bed', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')

            ##updating 111623
            ##check the loc_type
            final_type = 'Allcelltype'

            TE_loc = col[0] + '_' + col[1] + '_' + col[2]

            if final_type in store_type_dic:
                store_type_dic[final_type][TE_loc] = 1
            else:
                store_type_dic[final_type] = {}
                store_type_dic[final_type][TE_loc] = 1

    ##add the all information
    subfunction_summarize_te(store_type_dic, store_te_loc_nm_dic, store_te_prop_line_list)


    with open(opt_dir + '/opt_te_prop_loctype' + spe_prefix + '.txt', 'w+') as opt:
        for eachline in store_te_prop_line_list:
            opt.write(eachline + '\n')




def subfunction_make_acr_close_to_gene_cate (ipt_spe_gff, ipt_acr_fl,opt_dir,ipt_prefix):

    store_arc_line_list = []
    store_acr_celltype_dic = {}
    with open(ipt_acr_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            sumit = int((int(col[1]) + int(col[2])) / 2)
            final_line = col[0] + '\t' + str(sumit) + '\t' + str(sumit + 1) + '\t' + col[0] + '_' + col[1] + '_' + \
                         col[2]
            # final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[0] + '_' + col[1] + '_' + col[2]
            store_arc_line_list.append(final_line)

            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            mt = re.match('.+;(.+)',col[3])
            celltype = mt.group(1)
            if celltype == 'broadly_accessible':
                finalcelltype = 'Broad'
            else:
                finalcelltype = 'CT'

            store_acr_celltype_dic[acrnm] = finalcelltype


    with open(opt_dir + '/temp_' + ipt_prefix + '_acr_sumit_bed.txt', 'w+') as opt:
        for eachline in store_arc_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_prefix + '_acr_sumit_bed.txt > ' + \
          opt_dir + '/temp_' + ipt_prefix + '_acr_sumit_bed_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##generate a gene bed
    store_gene_bed_line_list = []
    with open(ipt_spe_gff, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if not eachline.startswith('#'):

                if not col[0].startswith('ChrMt') and not col[0].startswith('ChrPt'):

                    if 'scaffold' not in col[0]:
                        if 'Gm' in col[0]:
                            if col[0].startswith('Gm'):
                                if 'gene' in col[2]:
                                    final_line = col[0] + '\t' + col[3] + '\t' + col[4] + '\t' + col[0] + '_' + col[
                                        3] + '_' + col[4]
                                    store_gene_bed_line_list.append(final_line)
                        else:
                            if 'gene' in col[2]:
                                final_line = col[0] + '\t' + col[3] + '\t' + col[4] + '\t' + col[0] + '_' + col[
                                    3] + '_' + \
                                             col[4]
                                store_gene_bed_line_list.append(final_line)

    with open(opt_dir + '/temp_' + ipt_prefix + '_gene_bed.txt', 'w+') as opt:
        for eachline in store_gene_bed_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_prefix + '_gene_bed.txt > ' + opt_dir + '/temp_' + ipt_prefix + '_gene_bed_sorted.bed'
    subprocess.call(cmd, shell=True)

    print('begin to plot closest')
    ##use the bedtools
    cmd = 'bedtools closest -a ' + opt_dir + '/temp_' + ipt_prefix + '_acr_sumit_bed_sorted.txt'+ \
          ' -b ' + opt_dir + '/temp_' + ipt_prefix + '_gene_bed_sorted.bed' + \
          ' -d > ' + opt_dir + '/temp_' + ipt_prefix + '_acr_sumit_gene_closest.txt'
    subprocess.call(cmd, shell=True)

    ##calculate the distance
    store_final_line_list = []
    with open(opt_dir + '/temp_' + ipt_prefix + '_acr_sumit_gene_closest.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            final_line = col[3] + '\t' + col[-1]
            store_final_line_list.append(final_line)

    with open(opt_dir + '/opt_' + ipt_prefix + '_acr_sumit_distance_to_close_gene.txt',
              'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##decide the cate
    store_acr_nearbyGeneCate_distance_dic = {}
    with open(opt_dir + '/temp_' + ipt_prefix + '_acr_sumit_gene_closest.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[3]
            dist = col[-1]

            if int(dist) <= 0:
                cate = 'Genic'
            else:
                if int(dist) > 2000:
                    cate = 'Distal'
                else:
                    cate = 'Proximal'

            store_acr_nearbyGeneCate_distance_dic[acrnm] = {'cate':cate,'dist':str(dist)}

            #store_acr_nearbyGeneCate_dic[acrnm] = cate

    with open(opt_dir + '/opt_' + ipt_prefix + '_acr_toGeneCate.txt', 'w+') as opt:
        for eachline in store_acr_nearbyGeneCate_distance_dic:
            celltype = store_acr_celltype_dic[eachline]
            opt.write(eachline + '\t' + store_acr_nearbyGeneCate_distance_dic[eachline]['cate'] + '\t' + store_acr_nearbyGeneCate_distance_dic[eachline]['dist'] + '\t' + celltype + '\n')
































