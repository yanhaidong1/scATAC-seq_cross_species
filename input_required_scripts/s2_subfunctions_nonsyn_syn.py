#!/usr/bin/env python

##updating 011624 we will check the enrichment of TE
##updating 011024 we will check the nonsyn TE enrichment and compared to the syn TE enrichment
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

#################
##updating 010323
##we will add more analysis for the motif enrichment analysis
##what motifs enriched on the nonsyn compared to the syn across all species
def subfunction_store_motif_ACR(ipt_intersect_fl):
    store_motif_ACR_dic = {}
    with open(ipt_intersect_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acr = col[0] + '_' + col[1] + '_' + col[2]
            mt = re.match('(.+)_.+', col[6])
            motif = mt.group(1)

            if motif in store_motif_ACR_dic:
                store_motif_ACR_dic[motif][acr] = 1
            else:
                store_motif_ACR_dic[motif] = {}
                store_motif_ACR_dic[motif][acr] = 1

    return (store_motif_ACR_dic)


def multi_run(args):
    return subfunction_enrichment(*args)

def subfunction_enrichment (store_non_syn_spe1_celltype_loc_dic, store_syn_spe1_celltype_loc_dic, target_celltype,ipt_fimo_motif_prediction_fl, opt_dir,
                            spe1_prefix,spe2_prefix):

    store_cell_type_syn_acr_dic = store_syn_spe1_celltype_loc_dic[target_celltype]
    store_cell_type_non_syn_acr_dic = store_non_syn_spe1_celltype_loc_dic[target_celltype]


    with open(opt_dir + '/temp_syn_acr.txt', 'w+') as opt:
        for eachline in store_cell_type_syn_acr_dic:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_syn_acr.txt > ' + \
          opt_dir + '/temp_syn_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    with open(opt_dir + '/temp_non_syn_acr.txt', 'w+') as opt:
        for eachline in store_cell_type_non_syn_acr_dic:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_non_syn_acr.txt > ' + \
          opt_dir + '/temp_non_syn_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##intersect with the fimo motif file
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_syn_acr_sorted.txt ' + \
          ' -b ' + ipt_fimo_motif_prediction_fl + ' > ' + opt_dir + '/temp_syn_acr_motif_intersect.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_non_syn_acr_sorted.txt ' + \
          ' -b ' + ipt_fimo_motif_prediction_fl + ' > ' + opt_dir + '/temp_non_syn_acr_motif_intersect.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    store_motif_ACR_syn_dic = subfunction_store_motif_ACR(opt_dir + '/temp_syn_acr_motif_intersect.txt')
    store_motif_ACR_nonsyn_dic = subfunction_store_motif_ACR(opt_dir + '/temp_non_syn_acr_motif_intersect.txt')

    ##conduct the enrichment
    store_final_line_list = []
    for eachmotif in store_motif_ACR_nonsyn_dic:

        target_motif_nonsyn_num = len(list(store_motif_ACR_nonsyn_dic[eachmotif].keys()))

        target_motif_syn_num = 0
        if eachmotif in store_motif_ACR_syn_dic:
            target_motif_syn_num = len(list(store_motif_ACR_syn_dic[eachmotif].keys()))

        not_target_motif_nonsyn_acr_dic = {}
        for eachmotif2 in store_motif_ACR_nonsyn_dic:
            if eachmotif != eachmotif2:
                for eachacr in store_motif_ACR_nonsyn_dic[eachmotif2]:
                    not_target_motif_nonsyn_acr_dic[eachacr] = 1



        not_target_motif_nonsyn_acr_num = len(list(not_target_motif_nonsyn_acr_dic.keys()))

        non_target_motif_syn_acr_dic = {}
        for eachmotif2 in store_motif_ACR_syn_dic:
            if eachmotif != eachmotif2:
                for eachacr in store_motif_ACR_syn_dic[eachmotif2]:
                    non_target_motif_syn_acr_dic[eachacr] = 1


        non_target_motif_syn_acr_num = len(list(non_target_motif_syn_acr_dic.keys()))

        oddsratio, pvalue = stats.fisher_exact(
            [[target_motif_nonsyn_num, target_motif_syn_num],
             [not_target_motif_nonsyn_acr_num, non_target_motif_syn_acr_num]],
            alternative='greater')

        final_line = eachmotif + '\t' + str(pvalue) + '\t' + str(target_motif_nonsyn_num) + '\t' + str(
            target_motif_syn_num) + \
                     '\t' + str(not_target_motif_nonsyn_acr_num) + '\t' + str(non_target_motif_syn_acr_num)
        store_final_line_list.append(final_line)

    with open(opt_dir + '/opt_fisher_enrich_non_syn_motif_' + spe1_prefix + '_' + spe2_prefix + '_' + target_celltype + '.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



def subfunction_fisher_motif_enrichment_analysis_for_syn_nonsyn (ipt_fimo_motif_prediction_fl, ipt_summary_fl,
                                                                spe1_prefix, spe2_prefix, opt_dir):


    store_non_syn_spe1_celltype_loc_dic = {}
    store_syn_spe1_celltype_loc_dic = {}

    store_syn_acr_dic = {}
    store_non_syn_acr_dic = {}

    count = 0
    all_celltype_dic = {}
    with open(ipt_summary_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:

                ##store the different cate of acr
                synornonsyn = col[11]

                spe1_acrloc = '\t'.join(col[1].split('_'))

                spe1_celltype = col[8]

                ##for the not shared
                if synornonsyn == 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            all_celltype_dic[eachspe1_celltype] = 1

                            if eachspe1_celltype in store_non_syn_spe1_celltype_loc_dic:
                                store_non_syn_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_non_syn_spe1_celltype_loc_dic[eachspe1_celltype] = {}
                                store_non_syn_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                    store_non_syn_acr_dic[spe1_acrloc] = 1


                else:

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_syn_spe1_celltype_loc_dic:
                                store_syn_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_syn_spe1_celltype_loc_dic[eachspe1_celltype] = {}
                                store_syn_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                    store_syn_acr_dic[spe1_acrloc] = 1


    the_core_need = len(list(all_celltype_dic.keys()))

    ##create dir
    store_parallele_analysis_dir = opt_dir + '/store_parallele_analysis_dir'
    if not os.path.exists(store_parallele_analysis_dir):
        os.makedirs(store_parallele_analysis_dir)

    for x in range(0, int(the_core_need)):
        dir_code = x + 1
        temp_output_dir = store_parallele_analysis_dir + '/temp_output_' + str(dir_code) + 'dir'
        if not os.path.exists(temp_output_dir):
            os.makedirs(temp_output_dir)

    all_celltype_list = list(all_celltype_dic.keys())

    temp_all_output_dir_list = glob.glob(store_parallele_analysis_dir + '/*')


    ##subfunction_enrichment (store_non_syn_spe1_celltype_loc_dic, store_syn_spe1_celltype_loc_dic, target_celltype,ipt_fimo_motif_prediction_fl, opt_dir,
    #                        spe1_prefix,spe2_prefix)

    pool = Pool(int(the_core_need))
    run_list = []
    for x in range(0, int(the_core_need)):
        each_func_argument = (store_non_syn_spe1_celltype_loc_dic,
                              store_syn_spe1_celltype_loc_dic,
                              all_celltype_list[x],
                              ipt_fimo_motif_prediction_fl,
                              temp_all_output_dir_list[x],
                              spe1_prefix,spe2_prefix)
        run_list.append(each_func_argument)
    pool.map(multi_run, run_list)




#################
##updating 010324
##we will use the binomial test
##updating 120723
def subfunction_check_motif_enrichment_per_file (ipt_target_celltype_acr_loc_fl_list,ipt_total_spe_acr_fl,
                                                 ipt_motif_fl,s2_s4_repeat_times,opt_dir):

    store_combine_final_line_list = []
    for eachtarget_fl in ipt_target_celltype_acr_loc_fl_list:

        print('print target fl is ' + eachtarget_fl)

        mt = re.match('.+/(.+)',eachtarget_fl)
        flnm = mt.group(1)

        print('print flnm is ' + flnm)


        mt = re.match('opt_(.+)__(.+)_sorted\.txt',flnm)
        ipt_cate = mt.group(1)
        ipt_celltype = mt.group(2)

        store_control_dir = opt_dir + '/store_' + ipt_cate + '__' + ipt_celltype + '_control_dir'
        if not os.path.exists(store_control_dir):
            os.makedirs(store_control_dir)

        target_acrloc_dic = {}
        with open (eachtarget_fl,'r') as ipt:
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
        for i in range(int(s2_s4_repeat_times)):

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
        with open (opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_real_acr.txt','w+') as opt:
            for eachacrloc in list(target_acrloc_dic.keys()):
                acr_line = '\t'.join(eachacrloc.split('_'))
                opt.write(acr_line + '\n')

        total_target_acrloc_num = len(list(target_acrloc_dic.keys()))

        cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_real_acr.txt > ' + \
              opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_real_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##intersect with the motif
        cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_real_acr_sorted.txt' + ' -b ' + ipt_motif_fl + ' > ' + \
              opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_intersect_real_acr_motif.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        store_motif_real_acrloc_dic = {}
        with open (opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_intersect_real_acr_motif.txt','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                acrloc = col[0] + '_' + col[1] + '_' + col[2]
                mt = re.match('(.+)_.+',col[6])
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
            subprocess.call(cmd,shell=True)

            store_motif_control_acrloc_dic = {}
            with open(opt_dir + '/temp_' + ipt_cate + '__' + ipt_celltype + '_intersect_control_acr_motif.txt', 'r') as ipt:
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
                prop = len(list(acrloc_dic.keys()))/total_target_acrloc_num

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

            final_line = ipt_cate + '\t' + ipt_celltype + '\t' +  eachmotif + '\t' + str(pval) + '\t' + str(real_acr_num) + '\t' + str(total_target_acrloc_num) + '\t' + str(avg_prop)
            store_combine_final_line_list.append(final_line)


    with open(opt_dir + '/opt_final_cate_motif_enrichment.txt', 'w+') as opt:
        for eachline in store_combine_final_line_list:
            opt.write(eachline + '\n')




def multi_run_bionmial_test(args):
    return subfunction_check_motif_enrichment_per_file(*args)


def subfunction_check_motif_enrichment_per_celltype_binomial (ipt_motif_fl,ipt_final_summary_fl,ipt_total_spe_acr_fl,ipt_opt_dir,
                                                     s2_s4_repeat_times,input_core_num):

    ##1. we will divide them into three category
    ##2. for each category, how many ACR per motif
    ##3. use the binomial test to check the enrichment

    store_non_syn_spe1_celltype_loc_dic = {}
    store_syn_spe1_celltype_loc_dic = {}

    count = 0
    all_celltype_dic = {}
    with open(ipt_final_summary_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:

                ##store the different cate of acr
                synornonsyn = col[11]

                spe1_acrloc = '\t'.join(col[1].split('_'))

                spe1_celltype = col[8]

                ##for the not shared
                if synornonsyn == 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            all_celltype_dic[eachspe1_celltype] = 1

                            if eachspe1_celltype in store_non_syn_spe1_celltype_loc_dic:
                                store_non_syn_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_non_syn_spe1_celltype_loc_dic[eachspe1_celltype] = {}
                                store_non_syn_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                else:

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_syn_spe1_celltype_loc_dic:
                                store_syn_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_syn_spe1_celltype_loc_dic[eachspe1_celltype] = {}
                                store_syn_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1



    ##Here we will save the multiple files into a dir
    store_all_target_acrloc_dir = ipt_opt_dir + '/store_all_target_acrloc_dir'
    if not os.path.exists(store_all_target_acrloc_dir):
        os.makedirs(store_all_target_acrloc_dir)

    two_category_list = ['nonsyn','syn']
    for eachcate in two_category_list:

        target_loc_dic = {}
        if eachcate == 'nonsyn':
            target_loc_dic = store_non_syn_spe1_celltype_loc_dic
        if eachcate == 'syn':
            target_loc_dic = store_syn_spe1_celltype_loc_dic


        for eachspe1_celltype in target_loc_dic:

            spe1_acrloc_dic = target_loc_dic[eachspe1_celltype]

            store_final_line_list = []
            for eachloc in spe1_acrloc_dic:
                loc_line = '\t'.join(eachloc.split(','))
                store_final_line_list.append(loc_line)

            with open (store_all_target_acrloc_dir + '/opt_' + eachcate + '__' + eachspe1_celltype + '.txt','w+') as opt:
                for eachline in store_final_line_list:
                    opt.write(eachline + '\n')

            cmd = 'sort -k1,1V -k2,2n ' + store_all_target_acrloc_dir + '/opt_' + eachcate + '__' + eachspe1_celltype + '.txt > ' + \
                  store_all_target_acrloc_dir + '/opt_' + eachcate + '__' + eachspe1_celltype + '_sorted.txt'
            print(cmd)
            subprocess.call(cmd,shell=True)


    all_target_fl_list = glob.glob(store_all_target_acrloc_dir + '/*_sorted.txt')


    ##create dir
    store_parallele_analysis_dir = ipt_opt_dir + '/store_parallele_analysis_dir'
    if not os.path.exists(store_parallele_analysis_dir):
        os.makedirs(store_parallele_analysis_dir)

    for x in range(0, int(input_core_num)):
        dir_code = x + 1
        temp_output_dir = store_parallele_analysis_dir + '/temp_output_' + str(dir_code) + 'dir'
        if not os.path.exists(temp_output_dir):
            os.makedirs(temp_output_dir)

    store_core_dic = {}
    core_count = -1
    array_split_list = np.array_split(all_target_fl_list, int(input_core_num))
    for eacharray in array_split_list:
        core_count += 1

        ##the array list contains cluster information
        array_list = []
        for eachitem in eacharray:
            array_list.append(eachitem)

        store_core_dic[str(core_count)] = array_list

    temp_all_output_dir_list = glob.glob(store_parallele_analysis_dir + '/*')


    ##subfunction_check_motif_enrichment_per_file (ipt_target_celltype_acr_loc_fl_list,ipt_total_spe_acr_fl,
    #                                             ipt_motif_fl,s2_s4_repeat_times,opt_dir)

    pool = Pool(int(input_core_num))
    run_list = []
    for x in range(0, int(input_core_num)):
        each_func_argument = (store_core_dic[str(x)],
                              ipt_total_spe_acr_fl,
                              ipt_motif_fl,
                              s2_s4_repeat_times,
                              temp_all_output_dir_list[x])
        run_list.append(each_func_argument)
    pool.map(multi_run_bionmial_test, run_list)



def subfunction_collect_enrich_binomial_results (ipt_opt_dir,opt_dir,ipt_spe1_prefix,ipt_spe2_prefix):

    store_parallele_analysis_dir = ipt_opt_dir + '/store_parallele_analysis_dir'

    all_temp_dir_list = glob.glob(store_parallele_analysis_dir + '/*')


    store_all_enrichment_fl_list = []
    for eachtempdir in all_temp_dir_list:
        opt_final_cate_motif_enrichment_fl = eachtempdir + '/opt_final_cate_motif_enrichment.txt'
        store_all_enrichment_fl_list.append(opt_final_cate_motif_enrichment_fl)

    cmd = 'cat ' + ' '.join(store_all_enrichment_fl_list) + ' > ' + opt_dir + '/opt_' + ipt_spe1_prefix + '_' + ipt_spe2_prefix + '_binomial_motif_enrich.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)




#################
##updating 011024 we will check if we could find the enrichment of TEs
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

                                ##updating 011224 only divide it to the auto and nonauto
                                if 'DNAauto' in new_te_all_fam_nm:
                                    new_super_fam_nm = 'DNAauto'
                                else:
                                    if 'DNAnona' in new_te_all_fam_nm:
                                        new_super_fam_nm = 'DNAnona'
                                    else:
                                        new_super_fam_nm = 'DNAothers_' + new_te_all_fam_nm


                                #new_super_fam_nm = 'nMITE'

                        else:
                            new_super_fam_nm = 'DNAothers_' + new_te_all_fam_nm

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


def subfunction_intersect_ACR_TEs (ipt_target_ACR_fl, ipt_TE_loc_fl, opt_dir,opt_prefix):

    ##intersect with the synACR
    cmd = 'bedtools intersect -wa -wb -a ' + ipt_TE_loc_fl + ' -b ' + ipt_target_ACR_fl + ' > ' + \
          opt_dir + '/temp_TE_intersect_syn_ACR_' + opt_prefix + '.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    store_te_loc_nm_dic = {}
    with open(opt_dir + '/temp_TE_intersect_syn_ACR_' + opt_prefix + '.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            TE_loc = col[0] + '_' + col[1] + '_' + col[2]
            TE_nm = col[4]
            store_te_loc_nm_dic[TE_loc] = TE_nm

    store_type_dic = {}
    store_loctype_loc = {}
    with open(opt_dir + '/temp_TE_intersect_syn_ACR_' + opt_prefix + '.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')

            ##updating 111623
            ##check the loc_type
            #loc_type = col[8]

            #mt = re.match('.+;(.+)', loc_type)
            #celltype = mt.group(1)
            #if celltype == 'broadly_accessible':
            #    final_type = 'Broad'
            #else:
            #    final_type = 'CT'

            TE_loc = col[0] + '_' + col[1] + '_' + col[2]

            if opt_prefix in store_type_dic:
                store_type_dic[opt_prefix][TE_loc] = 1
            else:
                store_type_dic[opt_prefix] = {}
                store_type_dic[opt_prefix][TE_loc] = 1

            te_acr_loc = col[5] + '_' + col[6] + '_' + col[7]
            store_loctype_loc[te_acr_loc] = opt_prefix

    ##Here we will store the CT and all at same time
    store_te_prop_line_list = []
    subfunction_summarize_te(store_type_dic, store_te_loc_nm_dic, store_te_prop_line_list)

    return (store_te_prop_line_list)


def subfunction_check_TE_all_fam_ACR_composition_syn_nonsyn (ipt_spe1_TE_fl,ipt_final_summary_fl,opt_dir):

    ##first check syn and nonsyn the TE composition

    ##check broad and CT ACR in syn and not in syn for the TE composition

    ##if there are some TE family enriched on the nonsyn other than the syn
    ##what kind of TE are enriched on the epidermal nonsyn other than the syn

    store_syn_acr_loc_dic = {}
    store_nonsyn_acr_loc_dic = {}
    store_syn_CT_acr_loc_dic = {}
    store_nonsyn_CT_acr_loc_dic = {}
    store_syn_broad_acr_loc_dic = {}
    store_nonsyn_broad_acr_loc_dic = {}


    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1

            if count != 1:

                acrloc = col[1]
                syntenic_update = col[11]
                acrcate = col[8]

                ##for the syntnic regions
                if syntenic_update != 'none':
                    store_syn_acr_loc_dic[acrloc] = 1

                    if acrcate == 'broadly_accessible':
                        store_syn_broad_acr_loc_dic[acrloc] = 1
                    else:
                        store_syn_CT_acr_loc_dic[acrloc] = 1

                else:
                    store_nonsyn_acr_loc_dic[acrloc] = 1

                    if acrcate == 'broadly_accessible':
                        store_nonsyn_broad_acr_loc_dic[acrloc] = 1
                    else:
                        store_nonsyn_CT_acr_loc_dic[acrloc] = 1

    store_all_cate_fl_dic = {}
    with open (opt_dir + '/temp_syn_acr.txt','w+') as opt:
        for eachline in store_syn_acr_loc_dic:
            mt = re.match('(.+)_(\d+)_(\d+)',eachline)
            locstr = mt.group(1) + '\t' + mt.group(2) + '\t' + mt.group(3)
            opt.write(locstr + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_syn_acr.txt > ' + opt_dir + '/temp_syn_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_all_cate_fl_dic['syn_all'] = opt_dir + '/temp_syn_acr_sorted.txt'


    with open (opt_dir + '/temp_non_syn_acr.txt','w+') as opt:
        for eachline in store_nonsyn_acr_loc_dic:
            mt = re.match('(.+)_(\d+)_(\d+)',eachline)
            locstr = mt.group(1) + '\t' + mt.group(2) + '\t' + mt.group(3)
            opt.write(locstr + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_non_syn_acr.txt > ' + opt_dir + '/temp_non_syn_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_all_cate_fl_dic['nonsyn_all'] = opt_dir + '/temp_non_syn_acr_sorted.txt'


    with open(opt_dir + '/temp_syn_CT_acr.txt', 'w+') as opt:
        for eachline in store_syn_CT_acr_loc_dic:
            mt = re.match('(.+)_(\d+)_(\d+)',eachline)
            locstr = mt.group(1) + '\t' + mt.group(2) + '\t' + mt.group(3)
            opt.write(locstr + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_syn_CT_acr.txt > ' + opt_dir + '/temp_syn_CT_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    store_all_cate_fl_dic['syn_CT'] = opt_dir + '/temp_syn_CT_acr_sorted.txt'


    with open(opt_dir + '/temp_syn_broad_acr.txt', 'w+') as opt:
        for eachline in store_syn_broad_acr_loc_dic:
            mt = re.match('(.+)_(\d+)_(\d+)',eachline)
            locstr = mt.group(1) + '\t' + mt.group(2) + '\t' + mt.group(3)
            opt.write(locstr + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_syn_broad_acr.txt > ' + opt_dir + '/temp_syn_broad_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    store_all_cate_fl_dic['syn_broad'] = opt_dir + '/temp_syn_broad_acr_sorted.txt'

    with open(opt_dir + '/temp_nonsyn_CT_acr.txt', 'w+') as opt:
        for eachline in store_nonsyn_CT_acr_loc_dic:
            mt = re.match('(.+)_(\d+)_(\d+)',eachline)
            locstr = mt.group(1) + '\t' + mt.group(2) + '\t' + mt.group(3)
            opt.write(locstr + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_nonsyn_CT_acr.txt > ' + opt_dir + '/temp_nonsyn_CT_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    store_all_cate_fl_dic['nonsyn_CT'] = opt_dir + '/temp_nonsyn_CT_acr_sorted.txt'

    with open(opt_dir + '/temp_nonsyn_broad_acr.txt', 'w+') as opt:
        for eachline in store_nonsyn_broad_acr_loc_dic:
            mt = re.match('(.+)_(\d+)_(\d+)',eachline)
            locstr = mt.group(1) + '\t' + mt.group(2) + '\t' + mt.group(3)
            opt.write(locstr + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_nonsyn_broad_acr.txt > ' + opt_dir + '/temp_nonsyn_broad_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    store_all_cate_fl_dic['nonsyn_broad'] = opt_dir + '/temp_nonsyn_broad_acr_sorted.txt'


    ##we will check the TEs
    store_final_line = []
    with open(ipt_spe1_TE_fl, 'r') as ipt:
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

    ##now we will get the different class of ACR
    store_final_line_list = []
    for eachcate in store_all_cate_fl_dic:

        opt_prefix = eachcate
        ipt_target_fl = store_all_cate_fl_dic[eachcate]

        store_te_prop_line_list = subfunction_intersect_ACR_TEs(ipt_target_fl,  opt_dir + '/temp_TE_sorted.bed', opt_dir, opt_prefix)

        for eachline in store_te_prop_line_list:

            final_line = opt_prefix + '\t' + eachline
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_final_TE_prop_in_eachcate.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

##updating 011424
##we will get the summary to cell type specific syn and nonsyn
def subfunction_check_TE_all_fam_ACR_composition_syn_nonsyn_per_celltype (ipt_spe1_TE_fl, ipt_final_summary_fl, opt_dir):
    ##first check syn and nonsyn the TE composition

    ##check broad and CT ACR in syn and not in syn for the TE composition

    ##if there are some TE family enriched on the nonsyn other than the syn
    ##what kind of TE are enriched on the epidermal nonsyn other than the syn

    store_syn_celltype_acr_loc_dic = {}
    store_nonsyn_celltype_acr_loc_dic = {}

    count = 0
    with open(ipt_final_summary_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1

            if count != 1:

                acrloc = col[1]
                syntenic_update = col[11]
                acrcate = col[8]

                ##for the syntnic regions
                if syntenic_update != 'none':

                    if acrcate != 'broadly_accessible':

                        celltype_list = acrcate.split(',')
                        for eachcelltype in celltype_list:

                            if eachcelltype in store_syn_celltype_acr_loc_dic:
                                store_syn_celltype_acr_loc_dic[eachcelltype][acrloc] = 1
                            else:
                                store_syn_celltype_acr_loc_dic[eachcelltype] = {}
                                store_syn_celltype_acr_loc_dic[eachcelltype][acrloc] = 1

                else:

                    if acrcate != 'broadly_accessible':

                        celltype_list = acrcate.split(',')
                        for eachcelltype in celltype_list:

                            if eachcelltype in store_nonsyn_celltype_acr_loc_dic:
                                store_nonsyn_celltype_acr_loc_dic[eachcelltype][acrloc] = 1
                            else:
                                store_nonsyn_celltype_acr_loc_dic[eachcelltype] = {}
                                store_nonsyn_celltype_acr_loc_dic[eachcelltype][acrloc] = 1


    two_cate_list = ['syn','nonsyn']
    store_all_cate_fl_dic = {}
    for eachcate in two_cate_list:

        if eachcate == 'syn':
            target_acr_loc_dic = store_syn_celltype_acr_loc_dic
        else:
            target_acr_loc_dic = store_nonsyn_celltype_acr_loc_dic

        for eachcelltype in target_acr_loc_dic:

            acrloc_dic = target_acr_loc_dic[eachcelltype]

            store_final_line_list = []
            for eachloc in acrloc_dic:

                mt = re.match('(.+)_(\d+)_(\d+)', eachloc)
                loc_str = mt.group(1) + '\t' + mt.group(2) + '\t' + mt.group(3)

                #loc_str = '\t'.join(eachloc.split('_'))
                store_final_line_list.append(loc_str)

            with open (opt_dir + '/temp_' + eachcelltype + '_' + eachcate + '_acr.txt','w+') as opt:
                for eachline in store_final_line_list:
                    opt.write(eachline + '\n')

            cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + eachcelltype + '_' + eachcate + '_acr.txt > ' + \
                  opt_dir + '/temp_' + eachcelltype + '_' + eachcate + '_acr_sorted.txt'
            print(cmd)
            subprocess.call(cmd,shell=True)

            if eachcate in store_all_cate_fl_dic:
                store_all_cate_fl_dic[eachcate][eachcelltype] = opt_dir + '/temp_' + eachcelltype + '_' + eachcate + '_acr_sorted.txt'
            else:
                store_all_cate_fl_dic[eachcate] = {}
                store_all_cate_fl_dic[eachcate][eachcelltype] = opt_dir + '/temp_' + eachcelltype + '_' + eachcate + '_acr_sorted.txt'

    ##we will check the TEs
    store_final_line = []
    with open(ipt_spe1_TE_fl, 'r') as ipt:
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

    ##now we will get the different class of ACR
    store_final_line_list = []
    for eachcate in store_all_cate_fl_dic:

        for eachCT in store_all_cate_fl_dic[eachcate]:

            opt_prefix = eachcate + '__' + eachCT
            ipt_target_fl = store_all_cate_fl_dic[eachcate][eachCT]

            store_te_prop_line_list = subfunction_intersect_ACR_TEs(ipt_target_fl, opt_dir + '/temp_TE_sorted.bed', opt_dir,
                                                                    opt_prefix)

            for eachline in store_te_prop_line_list:
                final_line = opt_prefix + '\t' + eachline
                store_final_line_list.append(final_line)

    with open(opt_dir + '/opt_final_TE_prop_in_eachcate_celltypeVer.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')





################
#conduct the TE enrichment per cell type between syn and nonsyn
def subfunction_fisher_enrichment (target_dic,target_prefix1,target_prefix2,opt_dir,final_prefix):

    nonsyn_TEnm_count_dic = target_dic[target_prefix1]
    syn_TEnm_count_dic = target_dic[target_prefix2]

    store_final_line_list = []
    for eachTEfam in nonsyn_TEnm_count_dic:

        target_TEfam_nonsyn_count = nonsyn_TEnm_count_dic[eachTEfam]

        if eachTEfam in syn_TEnm_count_dic:
            target_TEfam_syn_count = syn_TEnm_count_dic[eachTEfam]
        else:
            target_TEfam_syn_count = 0

        not_target_TEfam_nonsyn_count = 0
        for eachTEfam2 in nonsyn_TEnm_count_dic:
            if eachTEfam2 != eachTEfam:
                not_target_TEfam_nonsyn_count = not_target_TEfam_nonsyn_count + int(nonsyn_TEnm_count_dic[eachTEfam2])

        not_target_TEfam_syn_count = 0
        for eachTEfam2 in syn_TEnm_count_dic:
            if eachTEfam2 != eachTEfam:
                not_target_TEfam_syn_count = not_target_TEfam_syn_count + int(syn_TEnm_count_dic[eachTEfam2])

        oddsratio, pvalue = stats.fisher_exact(
            [[int(target_TEfam_nonsyn_count), int(target_TEfam_syn_count)],
             [int(not_target_TEfam_nonsyn_count), int(not_target_TEfam_syn_count)]], alternative='greater')

        final_line = eachTEfam + '\t' + str(pvalue) + '\t' + str(target_TEfam_nonsyn_count) + '\t' + str(
            target_TEfam_syn_count) + '\t' + str(not_target_TEfam_nonsyn_count) + '\t' + str(not_target_TEfam_syn_count)
        store_final_line_list.append(final_line)

    with open(opt_dir + '/opt_nonsyn_vs_syn_TE_enrich_fam_' + final_prefix + '.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


def subfunction_conduct_enrichment_from_prop_fl_fisher (ipt_count_prop_fl,opt_dir):

    store_cate_TE_fam_count_SynVsNonSyn_dic = {}
    store_cate_TE_fam_count_broad_SynVsNonSyn_dic = {}
    store_cate_TE_fam_count_CT_SynVsNonSyn_dic = {}

    with open(ipt_count_prop_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            TEnm = col[2]
            TEcount = col[3]

            #if 'LTR' in TEnm and 'non_LTR' not in TEnm:

            #    if 'LTR.Copia' in TEnm:
            #        newTEnm = 'LTR.Copia'
            #    else:
            #        if 'LTR.Gypsy' in TEnm:
            #            newTEnm = 'LTR.Gypsy'
            #        else:
            #            newTEnm = 'LTR.Others'

            #else:
            #    newTEnm = TEnm

            newTEnm = TEnm

            cate = col[0]

            if cate == 'syn_all' or cate == 'nonsyn_all':

                if cate in store_cate_TE_fam_count_SynVsNonSyn_dic:
                    store_cate_TE_fam_count_SynVsNonSyn_dic[cate][newTEnm] = TEcount

                else:
                    store_cate_TE_fam_count_SynVsNonSyn_dic[cate] = {}
                    store_cate_TE_fam_count_SynVsNonSyn_dic[cate][newTEnm] = TEcount

            if cate == 'syn_broad' or cate == 'nonsyn_broad':

                if cate in store_cate_TE_fam_count_broad_SynVsNonSyn_dic:
                    store_cate_TE_fam_count_broad_SynVsNonSyn_dic[cate][newTEnm] = TEcount

                else:
                    store_cate_TE_fam_count_broad_SynVsNonSyn_dic[cate] = {}
                    store_cate_TE_fam_count_broad_SynVsNonSyn_dic[cate][newTEnm] = TEcount

            if cate == 'syn_CT' or cate == 'nonsyn_CT':

                if cate in store_cate_TE_fam_count_CT_SynVsNonSyn_dic:
                    store_cate_TE_fam_count_CT_SynVsNonSyn_dic[cate][newTEnm] = TEcount

                else:
                    store_cate_TE_fam_count_CT_SynVsNonSyn_dic[cate] = {}
                    store_cate_TE_fam_count_CT_SynVsNonSyn_dic[cate][newTEnm] = TEcount


    subfunction_fisher_enrichment(store_cate_TE_fam_count_SynVsNonSyn_dic, 'nonsyn_all', 'syn_all', opt_dir, 'all')
    subfunction_fisher_enrichment(store_cate_TE_fam_count_broad_SynVsNonSyn_dic, 'nonsyn_broad', 'syn_broad', opt_dir, 'broad')
    subfunction_fisher_enrichment(store_cate_TE_fam_count_CT_SynVsNonSyn_dic, 'nonsyn_CT', 'syn_CT', opt_dir,'CT')


##udpating 011424
def subfunction_conduct_enrichment_from_prop_fl_fisher_per_celltype (ipt_count_prop_fl,opt_dir):


    store_eachcelltype_cate_TE_fam_count_dic = {}

    with open(ipt_count_prop_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            TEnm = col[2]
            TEcount = col[3]

            #if 'LTR' in TEnm and 'non_LTR' not in TEnm:

            #    if 'LTR.Copia' in TEnm:
            #        newTEnm = 'LTR.Copia'
            #    else:
            #        if 'LTR.Gypsy' in TEnm:
            #            newTEnm = 'LTR.Gypsy'
            #        else:
            #            newTEnm = 'LTR.Others'

            #else:
            #    newTEnm = TEnm

            newTEnm = TEnm

            cate = col[0]

            mt = re.match('(.+)__(.+)',cate)
            synornonsyn = mt.group(1)
            celltype = mt.group(2)

            if celltype in store_eachcelltype_cate_TE_fam_count_dic:

                if synornonsyn in store_eachcelltype_cate_TE_fam_count_dic[celltype]:

                    store_eachcelltype_cate_TE_fam_count_dic[celltype][synornonsyn][newTEnm] = TEcount

                else:

                    store_eachcelltype_cate_TE_fam_count_dic[celltype][synornonsyn] = {}
                    store_eachcelltype_cate_TE_fam_count_dic[celltype][synornonsyn][newTEnm] = TEcount

            else:
                store_eachcelltype_cate_TE_fam_count_dic[celltype] = {}
                store_eachcelltype_cate_TE_fam_count_dic[celltype][synornonsyn] = {}
                store_eachcelltype_cate_TE_fam_count_dic[celltype][synornonsyn][newTEnm] = TEcount


    store_final_line_list = []
    for eachcelltype in store_eachcelltype_cate_TE_fam_count_dic:

        final_store_dic = store_eachcelltype_cate_TE_fam_count_dic[eachcelltype]

        subfunction_fisher_enrichment(final_store_dic, 'nonsyn', 'syn', opt_dir, 'CTver__' + eachcelltype)
        opt_celltype_enrichment_fl = opt_dir + '/opt_nonsyn_vs_syn_TE_enrich_fam_' + 'CTver__' + eachcelltype + '.txt'

        with open (opt_celltype_enrichment_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                final_line = eachcelltype + '\t' + eachline
                store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_nonsyn_vs_syn_TE_enrich_fam_combineAllCT.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




##updating 011624
##conduct the binomial for the target TE overlapping motifs
def subfunction_conduct_motif_enrichment_target_CT_TEfam (ipt_target_CT_TEfam_str,opt_dir_store_intersect_dir,ipt_total_spe_acr_fl,ipt_motif_fl,opt_dir,
                                                          s2_s11_sub_repeat_times):

    ipt_target_CT_TEfam_list = ipt_target_CT_TEfam_str.split(',')

    ipt_target_celltype_acr_loc_fl_list = []
    for eachstr in ipt_target_CT_TEfam_list:

        mt = re.match('(.+):(.+)',eachstr)
        target_ct = mt.group(1)
        target_TEfam = mt.group(2)

        ipt_target_fl = opt_dir_store_intersect_dir + '/temp_TE_intersect_syn_ACR_nonsyn__' + target_ct + '.txt'

        store_final_line_dic = {}
        with open (ipt_target_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()

                if col[4] == target_TEfam:

                    acrloc = col[5] + '\t' + col[6] + '\t' + col[7]
                    store_final_line_dic[acrloc] = 1

        with open (opt_dir + '/temp_target_ACR_' + target_ct + '_' + target_TEfam + '.txt','w+') as opt:
            for eachline in store_final_line_dic:
                opt.write(eachline + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_target_ACR_' + target_ct + '_' + target_TEfam + '.txt > ' + \
              opt_dir + '/opt_' + target_TEfam + '__' + target_ct + '_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ipt_target_celltype_acr_loc_fl_list.append( opt_dir + '/opt_' + target_TEfam + '__' + target_ct + '_sorted.txt')

        ##mt = re.match('opt_(.+)__(.+)_sorted\.txt',flnm)

    ##conduct the enrichment test



    subfunction_check_motif_enrichment_per_file (ipt_target_celltype_acr_loc_fl_list,ipt_total_spe_acr_fl,
                                                ipt_motif_fl,s2_s11_sub_repeat_times,opt_dir)



##updating 011724
##we will perform the DNA methylation surrounding the target ACR-LTR and LTR not overlapped with
def subfunction_conduct_DNA_methylation_target_CT_TEfam (input_all_three_methy_context_dir,ipt_target_CT_TEfam_str,opt_dir_store_intersect_dir,
                                                         ipt_spe1_TE_fl,ipt_total_spe_acr_fl,input_core_num,opt_dir):

    ipt_target_CT_TEfam_list = ipt_target_CT_TEfam_str.split(',')

    for eachstr in ipt_target_CT_TEfam_list:

        mt = re.match('(.+):(.+)', eachstr)
        target_ct = mt.group(1)
        target_TEfam = mt.group(2)

        ipt_target_fl = opt_dir_store_intersect_dir + '/temp_TE_intersect_syn_ACR_nonsyn__' + target_ct + '.txt'


        store_final_line_ACR_dic = {}
        store_final_line_TE_dic = {}
        with open(ipt_target_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()

                if col[4] == target_TEfam:
                    acrloc = col[5] + '\t' + col[6] + '\t' + col[7]
                    store_final_line_ACR_dic[acrloc] = 1

                    TEloc = col[0] + '\t' + col[1] + '\t' + col[2]
                    store_final_line_TE_dic[TEloc] = 1


        with open (opt_dir + '/temp_target_ACR_' + target_ct + '_' + target_TEfam + '.txt', 'w+') as opt:
            for eachline in store_final_line_ACR_dic:
                opt.write(eachline + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_target_ACR_' + target_ct + '_' + target_TEfam + '.txt > ' + \
              opt_dir + '/temp_target_ACR_' + target_ct + '_' + target_TEfam + '_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        with open (opt_dir + '/temp_target_TE_' + target_ct + '_' + target_TEfam + '.txt', 'w+') as opt:
            for eachline in store_final_line_TE_dic:
                opt.write(eachline + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_target_TE_' + target_ct + '_' + target_TEfam + '.txt > ' + \
              opt_dir + '/temp_target_TE_' + target_ct + '_' + target_TEfam + '_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)


        ##intersect with the synACR


        cmd = 'bedtools intersect -wa -wb -a ' + opt_dir_store_intersect_dir + '/temp_TE_sorted.bed' + ' -b ' + ipt_total_spe_acr_fl + ' > ' + \
              opt_dir + '/temp_TE_intersect_all_ACR.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)

        store_TE_intersect_ACR_dic = {}
        with open ( opt_dir + '/temp_TE_intersect_all_ACR.txt','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                TEloc = col[0] + '\t' + col[1] + '\t' + col[2]
                store_TE_intersect_ACR_dic[TEloc] = 1

        store_final_line_list = []
        with open (opt_dir_store_intersect_dir + '/temp_TE_sorted.bed', 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                TEloc = col[0] + '\t' + col[1] + '\t' + col[2]
                TEfam = col[4]

                if TEfam == target_TEfam:

                    if TEloc not in store_TE_intersect_ACR_dic:
                        store_final_line_list.append(TEloc)

        with open (opt_dir + '/temp_TE_' + target_TEfam + '_notOverlapACR.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_TE_' + target_TEfam + '_notOverlapACR.txt > ' + \
              opt_dir + '/temp_TE_' + target_TEfam + '_notOverlapACR_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)


        ##now we will do the plotting
        all_bw_fl_list = glob.glob(input_all_three_methy_context_dir + '/*bw*')
        for eachbwfl in all_bw_fl_list:
            mt = re.match('.+/(.+)',eachbwfl)
            flnm = mt.group(1)

            mt = re.match('leafmethy_v3\.bw\.(.+)',flnm)
            context_type = mt.group(1)

            ##plot meta plot for total, endosperm and not endosperm acr methyaltion comparision
            ##we will use the sorted
            cmd = 'computeMatrix reference-point' + \
                  ' --referencePoint center' + \
                  ' -a 2000 -b 2000' + \
                  ' -S ' + eachbwfl + ' ' + \
                  ' -R ' + opt_dir + '/temp_target_ACR_' + target_ct + '_' + target_TEfam + '_sorted.txt' + ' ' + \
                           opt_dir + '/temp_target_TE_' + target_ct + '_' + target_TEfam + '_sorted.txt' + ' ' + \
                           opt_dir + '/temp_TE_' + target_TEfam + '_notOverlapACR_sorted.txt' + ' ' + \
                  ' -p ' + input_core_num + \
                  ' --sortRegions descend' + \
                  ' -out ' + opt_dir + '/temp_combine_' + context_type + '_target_TE_ACR.tab.gz'
            print(cmd)
            subprocess.call(cmd, shell=True)


            ##plot heatmap for not endosperm acr
            method_list = ['nearest']

            for eachmethod in method_list:
                cmd = 'plotHeatmap' + \
                      ' -m ' +  opt_dir + '/temp_combine_' + context_type + '_target_TE_ACR.tab.gz' + \
                      ' -out ' + opt_dir + '/opt_' + eachmethod + '_' + context_type + '_target_TE_ACR.png' + \
                      ' --colorList \"#fff9bf,#3b528b\" --missingDataColor \"#FFF6EB\" --interpolationMethod ' + eachmethod + \
                      ' --heatmapHeight 15 --zMin 0 --zMax 1 --sortRegions descend' + \
                      ' --plotTitle \'' + context_type + '\'' + \
                      ' --regionsLabel \'ACR or TE\''
                print(cmd)
                subprocess.call(cmd, shell=True)


            ##plot the heatmap with different method
            method_list = ['nearest']

            for eachmethod in method_list:
                cmd = 'plotProfile' + \
                      ' -m ' +opt_dir + '/temp_combine_' + context_type + '_target_TE_ACR.tab.gz' + \
                      ' -out ' + opt_dir + '/opt_' + eachmethod + '_combine_' + context_type + '_target_TE_ACR_meta.pdf' + \
                      ' --plotFileFormat pdf' + \
                      ' --plotTitle \'' + context_type + '\''

                print(cmd)
                subprocess.call(cmd, shell=True)










































