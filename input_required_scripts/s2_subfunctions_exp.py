#!/usr/bin/env python

##this script is to check the expression

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





##updating 010824
def subfunction_check_celltype_specific_exp_per_gene (store_gene_celltype_cpm_dic,ipt_cutoff_overlap_avg):

    store_gene_final_celltype_exp_str_dic = {}
    for eachgene in store_gene_celltype_cpm_dic:

        celltype_cpm_dic = store_gene_celltype_cpm_dic[eachgene]

        celltype_float_cpm_dic = {}
        all_float_cpm_list = []
        for eachcelltype in celltype_cpm_dic:
            celltype_float_cpm_dic[eachcelltype] = float(celltype_cpm_dic[eachcelltype])
            all_float_cpm_list.append(float(celltype_cpm_dic[eachcelltype]))

        avg_cpm = mean(all_float_cpm_list)
        overavg_cutoff = avg_cpm * float(ipt_cutoff_overlap_avg)

        store_celltype_overlap_cutoff_list = []
        for eachcelltype in celltype_float_cpm_dic:
            float_cpm = celltype_float_cpm_dic[eachcelltype]

            if float_cpm >= overavg_cutoff:
                store_celltype_overlap_cutoff_list.append(eachcelltype)

        if store_celltype_overlap_cutoff_list != []:
            final_celltype_specific_exp_str = ','.join(store_celltype_overlap_cutoff_list)
        else:
            final_celltype_specific_exp_str = 'none'

        store_gene_final_celltype_exp_str_dic[eachgene] = final_celltype_specific_exp_str

    return (store_gene_final_celltype_exp_str_dic)


def subfunction_check_celltype_specific_exp_per_gene_return_highest_version (store_gene_celltype_cpm_dic):

    store_gene_final_celltype_exp_str_dic = {}
    for eachgene in store_gene_celltype_cpm_dic:

        celltype_cpm_dic = store_gene_celltype_cpm_dic[eachgene]

        celltype_float_cpm_dic = {}
        #all_float_cpm_list = []
        store_celltype_cpm_str_list = []
        for eachcelltype in celltype_cpm_dic:
            celltype_float_cpm_dic[eachcelltype] = float(celltype_cpm_dic[eachcelltype])
            #all_float_cpm_list.append(float(celltype_cpm_dic[eachcelltype]))

            final_line = eachcelltype + ':' + str(float(celltype_cpm_dic[eachcelltype]))
            store_celltype_cpm_str_list.append(final_line)

        max_celltype = max(celltype_float_cpm_dic, key=celltype_float_cpm_dic.get)

        store_celltype_cpm_str = ','.join(store_celltype_cpm_str_list)

        store_gene_final_celltype_exp_str_dic[eachgene] = {'max_ct':max_celltype,'exp_str':store_celltype_cpm_str}

    return (store_gene_final_celltype_exp_str_dic)

def subfunction_check_ortho_non_syn_synSS_class_gene_nearby_ACRs (ipt_spe1_spe2_orth_fl,ipt_cluster_annot_fl,ipt_cpm_snRNAseq_spe1_fl, ipt_cpm_snRNAseq_spe2_fl,opt_dir,
                                                                  s2_s3_target_celltype_addCpm_spe1_str,s2_s3_target_celltype_addCpm_spe2_str,s2_s10_cutoff_avg_fold,
                                                                  spe1_prefix,spe2_prefix):

    ##we will build a table to show
    ##rice orth \t cell-type-specific \t maize cell-type-specific
    store_cluster_celltype_dic = {}
    with open(ipt_cluster_annot_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            clusterID = col[0]
            annot = col[2]
            store_cluster_celltype_dic[clusterID] = annot



    if spe1_prefix == 'rice':

        store_gene_celltype_cpm_dic = {}
        store_order_clusterID_dic = {}
        store_all_gene_spe1_dic = {}
        count = 0
        with open(ipt_cpm_snRNAseq_spe1_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()

                count += 1
                if count == 1:

                    for i in range(len(col)):
                        store_order_clusterID_dic[i + 1] = store_cluster_celltype_dic[col[i]]

                else:
                    geneID = col[0]

                    store_celltype_cpmval_dic = {}
                    for i in range(1, len(col)):

                        celltype = store_order_clusterID_dic[i]
                        cpmval = col[i]

                        if celltype in s2_s3_target_celltype_addCpm_spe1_str.split(','):
                            store_celltype_cpmval_dic[celltype] = cpmval

                    store_gene_celltype_cpm_dic[geneID] = store_celltype_cpmval_dic

                    store_all_gene_spe1_dic[geneID] = 1


        store_gene_celltype_cpm_spe2_dic = {}
        store_order_celltype_spe2_dic = {}
        store_all_gene_spe2_dic = {}
        count = 0
        with open(ipt_cpm_snRNAseq_spe2_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1
                if count == 1:
                    for i in range(len(col)):
                        store_order_celltype_spe2_dic[i + 1] = col[i]

                else:
                    geneID = col[0]
                    store_celltype_cpmval_dic = {}
                    for i in range(1, len(col)):
                        celltype = store_order_celltype_spe2_dic[i]
                        cpmval = col[i]

                        if celltype in s2_s3_target_celltype_addCpm_spe2_str.split(','):
                            store_celltype_cpmval_dic[celltype] = cpmval

                    store_gene_celltype_cpm_spe2_dic[geneID] = store_celltype_cpmval_dic

                    store_all_gene_spe2_dic[geneID] = 1

    else:


        store_gene_celltype_cpm_dic = {}
        store_order_celltype_dic = {}
        store_all_gene_spe1_dic = {}
        count = 0
        with open(ipt_cpm_snRNAseq_spe1_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1
                if count == 1:
                    for i in range(len(col)):
                        store_order_celltype_dic[i + 1] = col[i]

                else:
                    geneID = col[0]
                    store_celltype_cpmval_dic = {}
                    for i in range(1, len(col)):
                        celltype = store_order_celltype_dic[i]
                        cpmval = col[i]

                        if celltype in s2_s3_target_celltype_addCpm_spe1_str.split(','):
                            store_celltype_cpmval_dic[celltype] = cpmval

                    store_gene_celltype_cpm_dic[geneID] = store_celltype_cpmval_dic

                    store_all_gene_spe1_dic[geneID] = 1


        store_gene_celltype_cpm_spe2_dic = {}
        store_order_clusterID_spe2_dic = {}
        store_all_gene_spe2_dic = {}
        count = 0
        with open(ipt_cpm_snRNAseq_spe2_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()

                count += 1
                if count == 1:

                    for i in range(len(col)):
                        store_order_clusterID_spe2_dic[i + 1] = store_cluster_celltype_dic[col[i]]

                else:
                    geneID = col[0]

                    store_celltype_cpmval_dic = {}
                    for i in range(1, len(col)):

                        celltype = store_order_clusterID_spe2_dic[i]
                        cpmval = col[i]

                        if celltype in s2_s3_target_celltype_addCpm_spe2_str.split(','):
                            store_celltype_cpmval_dic[celltype] = cpmval

                    store_gene_celltype_cpm_spe2_dic[geneID] = store_celltype_cpmval_dic

                    store_all_gene_spe2_dic[geneID] = 1

    ipt_cutoff_overlap_avg = s2_s10_cutoff_avg_fold

    #print(store_gene_celltype_cpm_dic)


    store_gene_final_celltype_exp_str_spe1_dic = subfunction_check_celltype_specific_exp_per_gene_return_highest_version (store_gene_celltype_cpm_dic)

    #store_gene_final_celltype_exp_str_spe1_dic = subfunction_check_celltype_specific_exp_per_gene (store_gene_celltype_cpm_dic,ipt_cutoff_overlap_avg)

    #store_gene_final_celltype_exp_str_spe2_dic = subfunction_check_celltype_specific_exp_per_gene (store_gene_celltype_cpm_spe2_dic,ipt_cutoff_overlap_avg)

    store_gene_final_celltype_exp_str_spe2_dic = subfunction_check_celltype_specific_exp_per_gene_return_highest_version (store_gene_celltype_cpm_spe2_dic)


    store_final_line_list = []
    count = 0
    with open (ipt_spe1_spe2_orth_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')

            count += 1
            if count != 1:

                if ',' not in col[1] and ',' not in col[2]:

                    if 'LOC' in col[1]:

                        mt = re.match('(.+)\.\d+',col[1])
                        spe1_gene = mt.group(1)

                        mt = re.match('(.+)_.+',col[2])
                        spe2_gene = mt.group(1)

                    else:

                        mt = re.match('(.+)_.+', col[1])
                        spe1_gene = mt.group(1)

                        mt = re.match('(.+)\.\d+', col[2])
                        spe2_gene = mt.group(1)



                    ##add the cell type str to the spe1 and spe2
                    if spe1_gene in store_all_gene_spe1_dic:

                        ##it means this gene could be detected
                        if spe1_gene in store_gene_final_celltype_exp_str_spe1_dic:
                            celltype_spe1_str = store_gene_final_celltype_exp_str_spe1_dic[spe1_gene]['max_ct']
                            celltype_cpm_spe1_str = store_gene_final_celltype_exp_str_spe1_dic[spe1_gene]['exp_str']
                        else:
                            celltype_spe1_str = 'error'
                            celltype_cpm_spe1_str = 'error'

                    else:
                        celltype_spe1_str = 'notDetectInCPMfile'
                        celltype_cpm_spe1_str = 'notDetectInCPMfile'

                    if spe2_gene in store_all_gene_spe2_dic:

                        if spe2_gene in store_gene_final_celltype_exp_str_spe2_dic:
                            celltype_spe2_str = store_gene_final_celltype_exp_str_spe2_dic[spe2_gene]['max_ct']
                            celltype_cpm_spe2_str = store_gene_final_celltype_exp_str_spe2_dic[spe2_gene]['exp_str']

                        else:
                            celltype_spe2_str = 'broad'
                            celltype_cpm_spe2_str = 'error'

                    else:
                        celltype_spe2_str = 'notDetectInCPMfile'
                        celltype_cpm_spe2_str = 'notDetectInCPMfile'


                    final_line = spe1_gene + '\t' + celltype_spe1_str + '\t' + celltype_cpm_spe1_str + '\t' + spe2_gene + '\t' + celltype_spe2_str + '\t' + celltype_cpm_spe2_str
                    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_ortho_gene_celltype_specific_exp' + spe1_prefix + '_' + spe2_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


def subfunction_add_orth_celltype_to_final_summary_fl (ipt_final_summary_fl,ipt_ortho_gene_celltype_specific_fl,spe1_prefix,spe2_prefix,opt_dir):

    store_spe1_to_spe2_dic = {}
    with open (ipt_ortho_gene_celltype_specific_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            spe1_gene = col[0]
            store_spe1_to_spe2_dic[spe1_gene] = col[1] + '\t' + col[2] + '\t' + col[3] + '\t' + col[4] + '\t' + col[5]


    store_final_line_list = []
    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1
            if count != 1:

                spe1_gene = col[19]

                if spe1_gene in store_spe1_to_spe2_dic:
                    final_line = eachline + '\t' + store_spe1_to_spe2_dic[spe1_gene]
                else:
                    final_line = eachline + '\t' + 'noOrtho' + '\t' + 'noOrtho' + '\t' + 'noOrtho' + '\t' + 'noOrtho' + '\t' + 'noOrtho'

                store_final_line_list.append(final_line)

            else:

                final_line = eachline + '\t' + spe1_prefix + '_highestCT' + '\t' + spe1_prefix + '_all_CTexp' + '\t' + spe2_prefix + '_ortho' + '\t' + spe2_prefix + '_highestCT' + '\t' + \
                             spe2_prefix + '_all_CTexp'
                store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_add_synOrNsyn_orth_exp_CTstr.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




    spe1_organct_list = []
    spe2_organct_list = []
    count = 0
    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_add_synOrNsyn_orth_exp_CTstr.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1
            if count == 2:
                organct_exp_spe1_str = col[22]
                organct_exp_spe1_list = organct_exp_spe1_str.split(',')

                for eachorganctexp in organct_exp_spe1_list:
                    #print(eachorganctexp)

                    mt = re.match('(.+):.+',eachorganctexp)
                    organct = mt.group(1)
                    spe1_organct_list.append(spe1_prefix + '_' + organct)

                organct_exp_spe2_str = col[25]
                organct_exp_spe2_list = organct_exp_spe2_str.split(',')

                for eachorganctexp in organct_exp_spe2_list:

                    mt = re.match('(.+):.+', eachorganctexp)
                    organct = mt.group(1)
                    spe2_organct_list.append(spe2_prefix + '_' + organct)


    spe1_organct_str = '\t'.join(spe1_organct_list)
    spe2_organct_str = '\t'.join(spe2_organct_list)

    store_final_line_multiple_exp_list = []
    count = 0
    with open(opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_add_synOrNsyn_orth_exp_CTstr.txt',
            'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1

            if count != 1:

                organct_exp_spe1_str = col[22]
                organct_exp_spe1_list = organct_exp_spe1_str.split(',')

                first_21_list = []
                for i in range(22):
                    first_21_list.append(col[i])

                if len(organct_exp_spe1_list) > 1:

                    spe1_exp_list = []
                    for eachorganctexp in organct_exp_spe1_list:
                        mt = re.match('.+:(.+)', eachorganctexp)
                        exp = mt.group(1)
                        spe1_exp_list.append(exp)

                    spe1_exp_str = '\t'.join(spe1_exp_list)

                else:
                    spe1_exp_str = 'none' + '\t' + 'none' + '\t' + 'none'


                organct_exp_spe2_str = col[25]
                organct_exp_spe2_list = organct_exp_spe2_str.split(',')

                if len(organct_exp_spe2_list) > 1:

                    spe2_exp_list = []
                    for eachorganctexp in organct_exp_spe2_list:
                        mt = re.match('.+:(.+)',eachorganctexp)
                        exp = mt.group(1)
                        spe2_exp_list.append(exp)

                    spe2_exp_str = '\t'.join(spe2_exp_list)

                else:
                    spe2_exp_str = 'none' + '\t' + 'none' + '\t' + 'none'

                final_line = '\t'.join(first_21_list) + '\t' + spe1_exp_str + '\t' + col[23] + '\t' + col[24] + '\t' + spe2_exp_str
                store_final_line_multiple_exp_list.append(final_line)

            else:
                first_21_list = []
                for i in range(22):
                    first_21_list.append(col[i])

                final_line = '\t'.join(first_21_list) + '\t' + spe1_organct_str + '\t' + col[23] + '\t' + col[24] + '\t' + spe2_organct_str
                store_final_line_multiple_exp_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_add_synOrNsyn_orth_exp_CTstr_collapes.txt','w+') as opt:
        for eachline in store_final_line_multiple_exp_list:
            opt.write(eachline + '\n')



##updating 011024
##we will extract the target acr loc from the final summary file
def subfunction_prepare_target_acr_loc (ipt_final_summary_fl,opt_dir,s2_s10_target_ct_str,ipt_spe1_prefix):

    s2_s10_target_ct_list = s2_s10_target_ct_str.split(',')

    for eachct in s2_s10_target_ct_list:

        if ipt_spe1_prefix == 'rice':

            store_spe2_not_target_ct = []
            store_spe1_target_ct = []
            if eachct == 'seedling.Epidermis':
                store_spe2_not_target_ct = ['mesophyll','companion_cells_sieve_elements']
                store_spe1_target_ct = ['protoderm','epidermis']
            if eachct == 'seedling.CompanionCell':
                store_spe2_not_target_ct = ['mesophyll', 'protoderm']
                store_spe1_target_ct = ['companion_cell']
            if eachct == 'seedling.Mesophyll':
                store_spe2_not_target_ct = ['companion_cells_sieve_elements', 'protoderm']
                store_spe1_target_ct = ['mesophyll']

        else:

            store_spe2_not_target_ct = []
            store_spe1_target_ct = []
            if eachct == 'protoderm':
                store_spe2_not_target_ct = ['seedling.Mesophyll', 'seedling.CompanionCell']
                store_spe1_target_ct = ['protoderm', 'epidermis']
            if eachct == 'companion_cells_sieve_elements':
                store_spe2_not_target_ct = ['seedling.Mesophyll', 'seedling.Epidermis']
                store_spe1_target_ct = ['companion_cells_sieve_elements']
            if eachct == 'mesophyll':
                store_spe2_not_target_ct = ['seedling.CompanionCell', 'seedling.Epidermis']
                store_spe1_target_ct = ['mesophyll']




        ##combine the syn SS and non syn regions
        store_target_spe1_acr_dic = {}
        with open (ipt_final_summary_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split('\t')

                syn_region = col[11]
                spe2_region = col[4]
                spe1_highestCT = col[21]
                spe2_highestCT = col[26]
                spe1_celltype = col[8]

                spe1_acr_loc = col[1]

                ##in the syn region
                if syn_region != 'none':

                    ##col4 is maize region
                    ##syn SS
                    if spe2_region == 'none':

                        if spe1_highestCT == eachct:

                            if spe2_highestCT in store_spe2_not_target_ct:

                                check_if_keep = ''
                                for eachtargetct in store_spe1_target_ct:
                                    if eachtargetct in spe1_celltype:
                                        check_if_keep = 'yes'

                                if check_if_keep == 'yes':
                                    store_target_spe1_acr_dic[spe1_acr_loc] = 1

                else:
                    if spe2_region == 'none':

                        if spe1_highestCT == eachct:

                            if spe2_highestCT in store_spe2_not_target_ct:

                                check_if_keep = ''
                                for eachtargetct in store_spe1_target_ct:
                                    if eachtargetct in spe1_celltype:
                                        check_if_keep = 'yes'

                                if check_if_keep == 'yes':
                                    store_target_spe1_acr_dic[spe1_acr_loc] = 1

        with open (opt_dir + '/temp_' + eachct + '_targetACR.txt','w+') as opt:
            for eachline in store_target_spe1_acr_dic:
                opt.write('\t'.join(eachline.split('_')) + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + eachct + '_targetACR.txt > ' + \
              opt_dir + '/opt_targetACR__' + eachct + '_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)



def subfunction_check_motif_enrichment(ipt_target_celltype_acr_loc_fl_list, ipt_total_spe_acr_fl,
                                                ipt_motif_fl, s2_s10_repeat_times, opt_dir):

    store_combine_final_line_list = []
    for eachtarget_fl in ipt_target_celltype_acr_loc_fl_list:

        print('print target fl is ' + eachtarget_fl)

        mt = re.match('.+/(.+)', eachtarget_fl)
        flnm = mt.group(1)

        print('print flnm is ' + flnm)

        mt = re.match('opt_(.+)__(.+)_sorted\.txt', flnm)
        ipt_cate = mt.group(1)
        ipt_celltype = mt.group(2)

        store_control_dir = opt_dir + '/store_' + ipt_cate + '__' + ipt_celltype + '_control_dir'
        if not os.path.exists(store_control_dir):
            os.makedirs(store_control_dir)

        target_acrloc_dic = {}
        with open(eachtarget_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                target_acrloc_dic['_'.join(col)] = 1

        store_all_acrloc_dic = {}
        with open(ipt_total_spe_acr_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                acrloc = col[0] + '_' + col[1] + '_' + col[2]
                store_all_acrloc_dic[acrloc] = 1

        ##build the control
        for i in range(int(s2_s10_repeat_times)):

            target_acrloc_len = len(list(target_acrloc_dic.keys()))
            store_all_acrloc_list = list(store_all_acrloc_dic.keys())

            rd_acr_list = random.sample(store_all_acrloc_list, target_acrloc_len)

            store_final_line_list = []
            for eachacr in rd_acr_list:
                acr_line = '\t'.join(eachacr.split('_'))
                store_final_line_list.append(acr_line)

            with open(store_control_dir + '/c' + str(i) + '.txt', 'w+') as opt:
                for eachline in store_final_line_list:
                    opt.write(eachline + '\n')

            cmd = 'sort -k1,1V -k2,2n ' + store_control_dir + '/c' + str(i) + '.txt > ' + \
                  store_control_dir + '/c' + str(i) + '_sorted.txt'
            print(cmd)
            subprocess.call(cmd, shell=True)

        ##intersect ACR to the motif
        ##for the target
        with open(opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_real_acr.txt', 'w+') as opt:
            for eachacrloc in list(target_acrloc_dic.keys()):
                acr_line = '\t'.join(eachacrloc.split('_'))
                opt.write(acr_line + '\n')

        total_target_acrloc_num = len(list(target_acrloc_dic.keys()))

        cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_real_acr.txt > ' + \
              opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_real_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)

        ##intersect with the motif
        cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_real_acr_sorted.txt' + ' -b ' + ipt_motif_fl + ' > ' + \
              opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_intersect_real_acr_motif.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)

        store_motif_real_acrloc_dic = {}
        with open(opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_intersect_real_acr_motif.txt',
                  'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                acrloc = col[0] + '_' + col[1] + '_' + col[2]
                mt = re.match('(.+)_.+', col[6])
                motifnm = mt.group(1)

                if motifnm in store_motif_real_acrloc_dic:
                    store_motif_real_acrloc_dic[motifnm][acrloc] = 1
                else:
                    store_motif_real_acrloc_dic[motifnm] = {}
                    store_motif_real_acrloc_dic[motifnm][acrloc] = 1

        all_control_fl_list = glob.glob(store_control_dir + '/*_sorted.txt')

        store_motif_control_acrprop_list_dic = {}
        rep_time = 0
        for eachcontrol_fl in all_control_fl_list:
            rep_time += 1

            cmd = 'bedtools intersect -wa -wb -a ' + eachcontrol_fl + ' -b ' + ipt_motif_fl + ' > ' + \
                  opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_intersect_control_acr_motif.txt'
            print(cmd)
            subprocess.call(cmd, shell=True)

            store_motif_control_acrloc_dic = {}
            with open(opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_intersect_control_acr_motif.txt',
                      'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split()
                    acrloc = col[0] + '_' + col[1] + '_' + col[2]
                    mt = re.match('(.+)_.+', col[6])
                    motifnm = mt.group(1)

                    if motifnm in store_motif_control_acrloc_dic:
                        store_motif_control_acrloc_dic[motifnm][acrloc] = 1
                    else:
                        store_motif_control_acrloc_dic[motifnm] = {}
                        store_motif_control_acrloc_dic[motifnm][acrloc] = 1

            ##calculate the prop for each motif
            ##total_target_acrloc_num

            for eachmotif in store_motif_control_acrloc_dic:
                acrloc_dic = store_motif_control_acrloc_dic[eachmotif]
                prop = len(list(acrloc_dic.keys())) / total_target_acrloc_num

                if eachmotif in store_motif_control_acrprop_list_dic:
                    store_motif_control_acrprop_list_dic[eachmotif].append(prop)
                else:
                    store_motif_control_acrprop_list_dic[eachmotif] = []
                    store_motif_control_acrprop_list_dic[eachmotif].append(prop)

        ##conduct the binomial
        for eachmotif in store_motif_real_acrloc_dic:

            real_acr_num = len(list(store_motif_real_acrloc_dic[eachmotif].keys()))

            if eachmotif in store_motif_control_acrprop_list_dic:
                avg_prop = mean(store_motif_control_acrprop_list_dic[eachmotif])
            else:
                avg_prop = 0

            pval = stats.binom_test(real_acr_num, n=total_target_acrloc_num, p=avg_prop, alternative='greater')

            final_line = ipt_cate + '\t' + ipt_celltype + '\t' + eachmotif + '\t' + str(pval) + '\t' + str(
                real_acr_num) + '\t' + str(total_target_acrloc_num) + '\t' + str(avg_prop)
            store_combine_final_line_list.append(final_line)

    with open(opt_dir + '/opt_final_cate_motif_enrichment.txt', 'w+') as opt:
        for eachline in store_combine_final_line_list:
            opt.write(eachline + '\n')


























