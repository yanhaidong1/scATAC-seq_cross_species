#!/usr/bin/env python


##updating 012424 we will use the restricted ACRs to have counting
##updating 010524 we will check the non-syntenic region for each cell type
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


##updating 120723
##build the CNS intersection for the lineage specific ACRs
def subfunction_check_CNS_overlap_density (ipt_spe1_cns_gff_fl,ipt_final_summary_fl,opt_dir):

    store_all_target_loc_fl_dir = opt_dir + '/store_all_target_loc_fl_dir'
    if not os.path.exists(store_all_target_loc_fl_dir):
        os.makedirs(store_all_target_loc_fl_dir)

    store_intersect_cns_fl_dir = opt_dir + '/store_intersect_cns_fl_dir'
    if not os.path.exists(store_intersect_cns_fl_dir):
        os.makedirs(store_intersect_cns_fl_dir)

    ##check the number of col in the gff file
    number_col = 0
    with open(ipt_spe1_cns_gff_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split('\t')
                number_col = len(col)

    if number_col > 4:

        store_final_line_list = []
        with open (ipt_spe1_cns_gff_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                if not eachline.startswith('#'):
                    col = eachline.strip().split()
                    cns_line = col[0] + '\t' + col[3] + '\t' + col[4]
                    store_final_line_list.append(cns_line)

    else:
        store_final_line_list = []
        with open (ipt_spe1_cns_gff_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                final_line = col[0] + '\t' + col[1] + '\t' + col[2]
                store_final_line_list.append(final_line)

    with open (opt_dir + '/temp_cns.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_cns.txt > ' + \
          opt_dir + '/opt_cns_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##we will consider several cate
    store_shared_A_spe1_broadorCT_loc_dic = {}
    store_shared_I_spe1_broadorCT_loc_dic = {}
    store_not_shared_spe1_broadorCT_loc_dic = {}

    store_syntenic_spe1_broadorCT_loc_dic = {}
    store_nonsyntenic_spe1_broadorCT_loc_dic = {}

    store_shared_I_spe1_CTcelltype_loc_dic = {}
    store_not_shared_spe1_CTcelltype_loc_dic = {}
    store_shared_A_spe1_CTcelltype_loc_dic = {}

    store_syntenic_spe1_CTcelltype_loc_dic = {}
    store_nonsyntenic_spe1_CTcelltype_loc_dic = {}

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
                spe2_acrloc = col[5]

                spe1_celltype = col[8]

                spe2_region = col[4]


                ##for the shared-A
                if spe1_acrloc != 'none' and spe2_acrloc != 'none':

                    if spe1_celltype == 'broadly_accessible':
                        broadorCT = 'Broad'
                    else:
                        broadorCT = 'CT'

                    if broadorCT in store_shared_A_spe1_broadorCT_loc_dic:
                        store_shared_A_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1
                    else:
                        store_shared_A_spe1_broadorCT_loc_dic[broadorCT] = {}
                        store_shared_A_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_shared_A_spe1_CTcelltype_loc_dic:
                                store_shared_A_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_shared_A_spe1_CTcelltype_loc_dic[eachspe1_celltype] = {}
                                store_shared_A_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1


                ##for the shared-I
                if spe1_acrID != 'none' and spe2_acrloc == 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_shared_I_spe1_CTcelltype_loc_dic:
                                store_shared_I_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_shared_I_spe1_CTcelltype_loc_dic[eachspe1_celltype] = {}
                                store_shared_I_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                        broadorCT = 'CT'

                    else:
                        broadorCT = 'Broad'

                    if broadorCT in store_shared_I_spe1_broadorCT_loc_dic:
                        store_shared_I_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1
                    else:
                        store_shared_I_spe1_broadorCT_loc_dic[broadorCT] = {}
                        store_shared_I_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1




                ##for the not shared
                if spe1_acrID == 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_not_shared_spe1_CTcelltype_loc_dic:
                                store_not_shared_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_not_shared_spe1_CTcelltype_loc_dic[eachspe1_celltype] = {}
                                store_not_shared_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                        broadorCT = 'CT'

                    else:
                        broadorCT = 'Broad'


                    if broadorCT in store_not_shared_spe1_broadorCT_loc_dic:
                        store_not_shared_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1
                    else:
                        store_not_shared_spe1_broadorCT_loc_dic[broadorCT] = {}
                        store_not_shared_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1


                ##updating 010524
                ##for the syntenic region
                ##store_syntenic_spe1_broadorCT_loc_dic
                ##store_syntenic_spe1_CTcelltype_loc_dic
                if spe2_region != 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_syntenic_spe1_CTcelltype_loc_dic:
                                store_syntenic_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_syntenic_spe1_CTcelltype_loc_dic[eachspe1_celltype] = {}
                                store_syntenic_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                        broadorCT = 'CT'

                    else:
                        broadorCT = 'Broad'

                    if broadorCT in store_syntenic_spe1_broadorCT_loc_dic:
                        store_syntenic_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1
                    else:
                        store_syntenic_spe1_broadorCT_loc_dic[broadorCT] = {}
                        store_syntenic_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1

                else:

                    ##store_nonsyntenic_spe1_broadorCT_loc_dic
                    ##store_nonsyntenic_spe1_CTcelltype_loc_dic

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_nonsyntenic_spe1_CTcelltype_loc_dic:
                                store_nonsyntenic_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_nonsyntenic_spe1_CTcelltype_loc_dic[eachspe1_celltype] = {}
                                store_nonsyntenic_spe1_CTcelltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                        broadorCT = 'CT'

                    else:
                        broadorCT = 'Broad'

                    if broadorCT in store_nonsyntenic_spe1_broadorCT_loc_dic:
                        store_nonsyntenic_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1
                    else:
                        store_nonsyntenic_spe1_broadorCT_loc_dic[broadorCT] = {}
                        store_nonsyntenic_spe1_broadorCT_loc_dic[broadorCT][spe1_acrloc] = 1







    ##save the loc information
    for eachcate in store_shared_A_spe1_broadorCT_loc_dic:

        with open (store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcate + '_acr.txt','w+') as opt:
            for eachloc in store_shared_A_spe1_broadorCT_loc_dic[eachcate]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcate + '_acr.txt > ' + \
              store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcate + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

    for eachcate in store_shared_I_spe1_broadorCT_loc_dic:

        with open(store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcate + '_acr.txt', 'w+') as opt:
            for eachloc in store_shared_I_spe1_broadorCT_loc_dic[eachcate]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcate + '_acr.txt > ' + \
              store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcate + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)

    for eachcate in store_not_shared_spe1_broadorCT_loc_dic:

        with open(store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcate + '_acr.txt', 'w+') as opt:
            for eachloc in store_not_shared_spe1_broadorCT_loc_dic[eachcate]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcate + '_acr.txt > ' + \
              store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcate + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)


    ##updating 010524
    for eachcate in store_syntenic_spe1_broadorCT_loc_dic:

        with open(store_all_target_loc_fl_dir + '/opt_syntenic_' + eachcate + '_acr.txt', 'w+') as opt:
            for eachloc in store_syntenic_spe1_broadorCT_loc_dic[eachcate]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_syntenic_' + eachcate + '_acr.txt > ' + \
              store_all_target_loc_fl_dir + '/opt_syntenic_' + eachcate + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)

    for eachcate in store_nonsyntenic_spe1_broadorCT_loc_dic:

        with open(store_all_target_loc_fl_dir + '/opt_nonsyntenic_' + eachcate + '_acr.txt', 'w+') as opt:
            for eachloc in store_nonsyntenic_spe1_broadorCT_loc_dic[eachcate]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_nonsyntenic_' + eachcate + '_acr.txt > ' + \
              store_all_target_loc_fl_dir + '/opt_nonsyntenic_' + eachcate + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)






    for eachcelltype in store_shared_I_spe1_CTcelltype_loc_dic:

        with open (store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcelltype + '_acr.txt','w+') as opt:
            for eachloc in store_shared_I_spe1_CTcelltype_loc_dic[eachcelltype]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcelltype + '_acr.txt > ' + \
              store_all_target_loc_fl_dir + '/opt_shared_I_' + eachcelltype + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd, shell=True)

    for eachcelltype in store_not_shared_spe1_CTcelltype_loc_dic:

        with open (store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcelltype + '_acr.txt','w+') as opt:
            for eachloc in store_not_shared_spe1_CTcelltype_loc_dic[eachcelltype]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcelltype + '_acr.txt' + \
              ' > ' + store_all_target_loc_fl_dir + '/opt_not_shared_' + eachcelltype + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

    for eachcelltype in store_shared_A_spe1_CTcelltype_loc_dic:

        with open (store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcelltype + '_acr.txt','w+') as opt:
            for eachloc in store_shared_A_spe1_CTcelltype_loc_dic[eachcelltype]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcelltype + '_acr.txt' + \
              ' > ' + store_all_target_loc_fl_dir + '/opt_shared_A_' + eachcelltype + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)


    ##updating 010524
    for eachcelltype in store_syntenic_spe1_CTcelltype_loc_dic:

        with open (store_all_target_loc_fl_dir + '/opt_syntenic_' + eachcelltype + '_acr.txt','w+') as opt:
            for eachloc in store_syntenic_spe1_CTcelltype_loc_dic[eachcelltype]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_syntenic_' + eachcelltype + '_acr.txt' + \
              ' > ' + store_all_target_loc_fl_dir + '/opt_syntenic_' + eachcelltype + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)


    for eachcelltype in store_nonsyntenic_spe1_CTcelltype_loc_dic:

        with open (store_all_target_loc_fl_dir + '/opt_nonsyntenic_' + eachcelltype + '_acr.txt','w+') as opt:
            for eachloc in store_nonsyntenic_spe1_CTcelltype_loc_dic[eachcelltype]:
                loc_line = '\t'.join(eachloc.split('_'))
                opt.write(loc_line + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_all_target_loc_fl_dir + '/opt_nonsyntenic_' + eachcelltype + '_acr.txt' + \
              ' > ' + store_all_target_loc_fl_dir + '/opt_nonsyntenic_' + eachcelltype + '_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)





    all_target_fl_list = glob.glob(store_all_target_loc_fl_dir + '/*sorted.txt')

    store_final_line_list = []
    for eachfl in all_target_fl_list:

        mt = re.match('.+/(.+)',eachfl)
        flnm = mt.group(1)

        mt = re.match('opt_(.+)_acr_sorted\.txt',flnm)
        catenm = mt.group(1)

        store_all_acr_loc_dic = {}
        with open (eachfl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                acrloc = col[0] + '_' + col[1] + '_' + col[2]
                store_all_acr_loc_dic[acrloc] = 1

        cmd = 'bedtools intersect -wa -wb -a ' + eachfl + \
              ' -b ' + opt_dir + '/opt_cns_sorted.txt > ' + \
              opt_dir + '/opt_intersect_' + catenm + '_acr_cns.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##calculate the num of cns per ACR
        store_acr_cns_num_dic = {}
        with open (opt_dir + '/opt_intersect_' + catenm + '_acr_cns.txt','r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                acrloc = col[0] + '_' + col[1] + '_' + col[2]
                cnsloc = col[3] + '_' + col[4] + '_' + col[5]

                if acrloc in store_acr_cns_num_dic:
                    store_acr_cns_num_dic[acrloc][cnsloc] = 1
                else:
                    store_acr_cns_num_dic[acrloc] = {}
                    store_acr_cns_num_dic[acrloc][cnsloc] = 1


        for eachacr in store_all_acr_loc_dic:

            if eachacr in store_acr_cns_num_dic:
                cnsloc_num = len(list(store_acr_cns_num_dic[eachacr].keys()))
            else:
                cnsloc_num = 0

            final_line = catenm + '\t' + eachacr + '\t' + str(cnsloc_num)
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_final_cate_acr_cns_num.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')





def subfunction_check_enrichment_for_celltype (ipt_final_cate_acr_cns_num_fl,opt_dir):


    store_sharedI_celltype_finalcatenum_dic = {}
    store_not_shared_celltype_finalcatenum_dic = {}
    store_sharedA_celltype_finalcatenum_dic = {}
    with open (ipt_final_cate_acr_cns_num_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            catenm = col[0]
            cnsloc_num = col[2]

            if cnsloc_num != '0':
                final_cate = 'CoverCNS'
            else:
                final_cate = 'NotCoverCNS'

            if '_Broad' not in catenm and '_CT' not in catenm:

                if '_unknown_cells' not in catenm:

                    ##focuse on the shared_I
                    if 'shared_I_' in catenm:

                        if catenm in store_sharedI_celltype_finalcatenum_dic:

                            if final_cate in store_sharedI_celltype_finalcatenum_dic[catenm]:
                                store_sharedI_celltype_finalcatenum_dic[catenm][final_cate] += 1
                            else:
                                store_sharedI_celltype_finalcatenum_dic[catenm][final_cate] = 1

                        else:
                            store_sharedI_celltype_finalcatenum_dic[catenm] = {}
                            store_sharedI_celltype_finalcatenum_dic[catenm][final_cate] = 1


                    if 'not_shared' in catenm:

                        if catenm in store_not_shared_celltype_finalcatenum_dic:

                            if final_cate in store_not_shared_celltype_finalcatenum_dic[catenm]:
                                store_not_shared_celltype_finalcatenum_dic[catenm][final_cate] += 1
                            else:
                                store_not_shared_celltype_finalcatenum_dic[catenm][final_cate] = 1

                        else:
                            store_not_shared_celltype_finalcatenum_dic[catenm] = {}
                            store_not_shared_celltype_finalcatenum_dic[catenm][final_cate] = 1

                    if 'shared_A_' in catenm:

                        if catenm in store_sharedA_celltype_finalcatenum_dic:

                            if final_cate in store_sharedA_celltype_finalcatenum_dic[catenm]:
                                store_sharedA_celltype_finalcatenum_dic[catenm][final_cate] += 1
                            else:
                                store_sharedA_celltype_finalcatenum_dic[catenm][final_cate] = 1

                        else:
                            store_sharedA_celltype_finalcatenum_dic[catenm] = {}
                            store_sharedA_celltype_finalcatenum_dic[catenm][final_cate] = 1






    ##for the sharedI
    store_fisher_test_final_line_list = []
    store_plot_prop_final_line_list = []
    for target_celltype in store_sharedI_celltype_finalcatenum_dic:

        target_celltype_CoverCNS_num = store_sharedI_celltype_finalcatenum_dic[target_celltype]['CoverCNS']
        target_celltype_NotCoverCNS_num = store_sharedI_celltype_finalcatenum_dic[target_celltype]['NotCoverCNS']

        other_celltype_coverCNS_total_num = 0
        other_celltype_NotcoverCNS_total_num = 0
        for other_celltype in store_sharedI_celltype_finalcatenum_dic:

            if other_celltype != target_celltype:

                other_celltype_CoverCNS_num = store_sharedI_celltype_finalcatenum_dic[other_celltype]['CoverCNS']
                other_celltype_NotCoverCNS_num = store_sharedI_celltype_finalcatenum_dic[other_celltype]['NotCoverCNS']

                other_celltype_coverCNS_total_num = other_celltype_coverCNS_total_num + other_celltype_CoverCNS_num
                other_celltype_NotcoverCNS_total_num = other_celltype_NotcoverCNS_total_num + other_celltype_NotCoverCNS_num


        ##use the fisher
        oddsratio, pvalue = stats.fisher_exact([[target_celltype_CoverCNS_num, target_celltype_NotCoverCNS_num],
                                                [other_celltype_coverCNS_total_num, other_celltype_NotcoverCNS_total_num]],
                                               alternative='greater')


        final_line = 'shared_I' + '\t' + target_celltype + '\t' + str(pvalue) + '\t' + str(target_celltype_CoverCNS_num) + '\t' + str(target_celltype_NotCoverCNS_num) + '\t' + \
                       str(other_celltype_coverCNS_total_num) + '\t' + str(other_celltype_NotcoverCNS_total_num)

        store_fisher_test_final_line_list.append(final_line)


        ##build the proprotion
        final_line = 'shared_I' + '\t' + target_celltype + '\t' + 'TargetCT' + '\t' + str(target_celltype_CoverCNS_num/ (target_celltype_CoverCNS_num + target_celltype_NotCoverCNS_num)) + '\n' + \
                     'shared_I' + '\t' + target_celltype + '\t' + 'NotTargetCT' + '\t' + str(other_celltype_coverCNS_total_num/(other_celltype_coverCNS_total_num + other_celltype_NotcoverCNS_total_num))
        store_plot_prop_final_line_list.append(final_line)

    ##for the not shared
    for target_celltype in store_not_shared_celltype_finalcatenum_dic:

        target_celltype_CoverCNS_num = store_not_shared_celltype_finalcatenum_dic[target_celltype]['CoverCNS']
        target_celltype_NotCoverCNS_num = store_not_shared_celltype_finalcatenum_dic[target_celltype]['NotCoverCNS']

        other_celltype_coverCNS_total_num = 0
        other_celltype_NotcoverCNS_total_num = 0
        for other_celltype in store_not_shared_celltype_finalcatenum_dic:

            if other_celltype != target_celltype:
                other_celltype_CoverCNS_num = store_not_shared_celltype_finalcatenum_dic[other_celltype]['CoverCNS']
                other_celltype_NotCoverCNS_num = store_not_shared_celltype_finalcatenum_dic[other_celltype]['NotCoverCNS']

                other_celltype_coverCNS_total_num = other_celltype_coverCNS_total_num + other_celltype_CoverCNS_num
                other_celltype_NotcoverCNS_total_num = other_celltype_NotcoverCNS_total_num + other_celltype_NotCoverCNS_num

        ##use the fisher
        oddsratio, pvalue = stats.fisher_exact([[target_celltype_CoverCNS_num, target_celltype_NotCoverCNS_num],
                                                [other_celltype_coverCNS_total_num,
                                                 other_celltype_NotcoverCNS_total_num]],
                                               alternative='greater')

        final_line = 'not_shared' + '\t' + target_celltype + '\t' + str(pvalue) + '\t' + str(
            target_celltype_CoverCNS_num) + '\t' + str(target_celltype_NotCoverCNS_num) + '\t' + \
                     str(other_celltype_coverCNS_total_num) + '\t' + str(other_celltype_NotcoverCNS_total_num)

        store_fisher_test_final_line_list.append(final_line)

        ##build the proprotion
        final_line = 'not_shared' + '\t' + target_celltype + '\t' + 'TargetCT' + '\t' + str(
            target_celltype_CoverCNS_num / (target_celltype_CoverCNS_num + target_celltype_NotCoverCNS_num)) + '\n' + \
                     'not_shared' + '\t' + target_celltype + '\t' + 'NotTargetCT' + '\t' + str(other_celltype_coverCNS_total_num / (
                    other_celltype_coverCNS_total_num + other_celltype_NotcoverCNS_total_num))
        store_plot_prop_final_line_list.append(final_line)


    ##for the shared-A
    for target_celltype in store_sharedA_celltype_finalcatenum_dic:

        target_celltype_CoverCNS_num = store_sharedA_celltype_finalcatenum_dic[target_celltype]['CoverCNS']
        target_celltype_NotCoverCNS_num = store_sharedA_celltype_finalcatenum_dic[target_celltype]['NotCoverCNS']

        other_celltype_coverCNS_total_num = 0
        other_celltype_NotcoverCNS_total_num = 0
        for other_celltype in store_sharedA_celltype_finalcatenum_dic:

            if other_celltype != target_celltype:
                other_celltype_CoverCNS_num = store_sharedA_celltype_finalcatenum_dic[other_celltype]['CoverCNS']
                other_celltype_NotCoverCNS_num = store_sharedA_celltype_finalcatenum_dic[other_celltype]['NotCoverCNS']

                other_celltype_coverCNS_total_num = other_celltype_coverCNS_total_num + other_celltype_CoverCNS_num
                other_celltype_NotcoverCNS_total_num = other_celltype_NotcoverCNS_total_num + other_celltype_NotCoverCNS_num

        ##use the fisher
        oddsratio, pvalue = stats.fisher_exact([[target_celltype_CoverCNS_num, target_celltype_NotCoverCNS_num],
                                                [other_celltype_coverCNS_total_num,
                                                 other_celltype_NotcoverCNS_total_num]],
                                               alternative='greater')

        final_line = 'shared_A' + '\t' + target_celltype + '\t' + str(pvalue) + '\t' + str(
            target_celltype_CoverCNS_num) + '\t' + str(target_celltype_NotCoverCNS_num) + '\t' + \
                     str(other_celltype_coverCNS_total_num) + '\t' + str(other_celltype_NotcoverCNS_total_num)

        store_fisher_test_final_line_list.append(final_line)

        ##build the proprotion
        final_line = 'shared_A' + '\t' + target_celltype + '\t' + 'TargetCT' + '\t' + str(
            target_celltype_CoverCNS_num / (target_celltype_CoverCNS_num + target_celltype_NotCoverCNS_num)) + '\n' + \
                     'shared_A' + '\t' + target_celltype + '\t' + 'NotTargetCT' + '\t' + str(
            other_celltype_coverCNS_total_num / (
                        other_celltype_coverCNS_total_num + other_celltype_NotcoverCNS_total_num))
        store_plot_prop_final_line_list.append(final_line)



    with open (opt_dir + '/opt_fisher_test_allcelltype_enrich_on_CNS.txt','w+') as opt:
        for eachline in store_fisher_test_final_line_list:
            opt.write(eachline + '\n')

    with open (opt_dir + '/opt_plot_R_prop_allcelltype_enrich_on_CNS.txt','w+') as opt:
        for eachline in store_plot_prop_final_line_list:
            opt.write(eachline + '\n')


    ##for the CNS check the shared-A
    ##what happens to these ACRs


##updating 010824
#def subfunction_add_CNS_back_to_summary_fl(ipt_spe1_final_summary_fl, ipt_spe1_acr_cns_num_fl):
#    store_acr_cns_dic = {}
#    with open(ipt_spe1_acr_cns_num_fl, 'r') as ipt:


##upadting 121923
def subfunction_check_varACRotherSpe_inRiceRegion_in_Atlas (ipt_final_summary_fl,ipt_rice_atlas_fl,ipt_spe2_cns_gff_fl,input_rice_atlas_acr_coverage_fl,
                                                            opt_dir,spe1_prefix,spe2_prefix,s2_s2_organct_coverage_cutoff):


    ##save the cell type num for the spe2
    ##we will check the H3K27m3 for each of organ one by one
    store_all_rice_atlas_peak_dic = {}
    with open (ipt_rice_atlas_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            store_all_rice_atlas_peak_dic[acrnm] = 1


    store_peak_organct_val_dic = {}
    store_total_organct_dic = {}
    with open(input_rice_atlas_acr_coverage_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            peaknm = col[1]
            organct = col[0]
            coverage_val = col[2]
            posneg_cate = col[3]
            if peaknm in store_peak_organct_val_dic:
                store_peak_organct_val_dic[peaknm][organct] = {'coverage_val': coverage_val,
                                                               'posneg_cate': posneg_cate}
            else:
                store_peak_organct_val_dic[peaknm] = {}
                store_peak_organct_val_dic[peaknm][organct] = {'coverage_val': coverage_val,
                                                               'posneg_cate': posneg_cate}

            store_total_organct_dic[organct] = 1

    store_final_line_list = []
    store_eachpeak_celltype_str_num_dic = {}
    for eachpeak in store_peak_organct_val_dic:

        ##we want to allow the peaks are in the final peaks we called
        if eachpeak in store_all_rice_atlas_peak_dic:

            store_celltype_list = []
            store_cover_val_list = []

            for eachorganct in store_peak_organct_val_dic[eachpeak]:
                coverage_val = store_peak_organct_val_dic[eachpeak][eachorganct]['coverage_val']
                posneg_cate = store_peak_organct_val_dic[eachpeak][eachorganct]['posneg_cate']

                if float(coverage_val) >= float(s2_s2_organct_coverage_cutoff):

                    ##updating 050423
                    if posneg_cate == 'Positive':
                        store_celltype_list.append(eachorganct)
                        store_cover_val_list.append(coverage_val)

            store_celltype_str = ','.join(store_celltype_list)
            store_cover_val_str = ','.join(store_cover_val_list)

            if store_celltype_str != '':
                celltypenum = len(store_celltype_list)
                final_line = 'CutoffCover:' + s2_s2_organct_coverage_cutoff + '\t' + eachpeak + '\t' + store_celltype_str + '\t' + store_cover_val_str + '\t' + str(
                    celltypenum)
                store_final_line_list.append(final_line)
                store_eachpeak_celltype_str_num_dic[eachpeak] = {'celltypestr': store_celltype_str,
                                                                 'celltypenum': str(celltypenum)}



    ##check the number of col in the gff file
    number_col = 0
    with open(ipt_spe2_cns_gff_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split('\t')
                number_col = len(col)

    if number_col > 4:

        store_final_line_list = []
        with open(ipt_spe2_cns_gff_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                if not eachline.startswith('#'):
                    col = eachline.strip().split()
                    cns_line = col[0] + '\t' + col[3] + '\t' + col[4]
                    store_final_line_list.append(cns_line)

    else:
        store_final_line_list = []
        with open(ipt_spe2_cns_gff_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                final_line = col[0] + '\t' + col[1] + '\t' + col[2]
                store_final_line_list.append(final_line)

    with open(opt_dir + '/temp_' + spe2_prefix + '_cns.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe2_prefix + '_cns.txt > ' + \
          opt_dir + '/opt_' + spe2_prefix + '_cns_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)




    store_blastedRegion_dic_in_rice_dic = {}
    store_spe1_acr_dic = {}
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
                spe2_acrloc = col[5]
                spe2_region = col[4]

                if spe2_region != 'none':
                    spe2_region_loc = '\t'.join(spe2_region.split('_'))
                    store_blastedRegion_dic_in_rice_dic[spe2_region_loc] = 1

                spe1_acr = '\t'.join(spe1_acrloc.split('_'))
                store_spe1_acr_dic[spe1_acr] = 1

    with open (opt_dir + '/temp_' + spe1_prefix + '_acr.txt','w+') as opt:
        for eachline in store_spe1_acr_dic:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt > ' + \
          opt_dir + '/temp_' + spe1_prefix + '_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)


    with open (opt_dir + '/temp_' + spe2_prefix+ '_region.txt','w+') as opt:
        for eachline in store_blastedRegion_dic_in_rice_dic:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe2_prefix + '_region.txt > ' + opt_dir + '/temp_' + spe2_prefix + '_region_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##intersect with the Atlas
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + spe2_prefix + '_region_sorted.txt' + \
          ' -b ' + ipt_rice_atlas_fl + ' > ' + opt_dir + '/temp_' + spe2_prefix + '_region_intersect_atlas_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_acrregion_atlas_acr_celltype_dic = {}
    with open (opt_dir + '/temp_' + spe2_prefix + '_region_intersect_atlas_acr.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            region_loc = col[0] + '_' + col[1] + '_' + col[2]
            rice_atlas_acr = col[3] + '_' + col[4] + '_' + col[5]
            celltype = col[7]

            if region_loc in store_acrregion_atlas_acr_celltype_dic:
                store_acrregion_atlas_acr_celltype_dic[region_loc][rice_atlas_acr] = celltype
            else:
                store_acrregion_atlas_acr_celltype_dic[region_loc] = {}
                store_acrregion_atlas_acr_celltype_dic[region_loc][rice_atlas_acr] = celltype


    ##opt_dir + '/opt_cns_sorted.txt'
    #intersect the cns with the spe1 ACR and spe2 region
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + spe1_prefix + '_acr_sorted.txt' + \
          ' -b ' + opt_dir + '/opt_cns_sorted.txt > ' + opt_dir + '/temp_' + spe1_prefix + '_intersect_acr_cns.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_spe1_acr_capture_cns_dic = {}
    with open (opt_dir + '/temp_' + spe1_prefix + '_intersect_acr_cns.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            store_spe1_acr_capture_cns_dic[acrloc] = 1


    ##intersect the spe2 cns with the sp2 region
    cmd = 'bedtools intersect -wa -wb -a ' +  opt_dir + '/temp_' + spe2_prefix + '_region_sorted.txt' + \
          ' -b ' + opt_dir + '/opt_' + spe2_prefix + '_cns_sorted.txt > ' + \
          opt_dir + '/temp_intersect_' + spe2_prefix + '_region_cns.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_spe2_region_contain_cns_dic = {}
    with open (opt_dir + '/temp_intersect_' + spe2_prefix + '_region_cns.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            regionloc = col[0] + '_' + col[1] + '_' + col[2]
            store_spe2_region_contain_cns_dic[regionloc] = 1


    store_final_line_list = []
    with open (ipt_final_summary_fl,'r') as ipt:
        count = 0
        with open(ipt_final_summary_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1
                if count != 1:
                    ##store the different cate of acr
                    spe2_region = col[4]
                    spe1_acr = col[1]

                    if spe1_acr in store_spe1_acr_capture_cns_dic:
                        spe1_cover_cns = 'yes'
                    else:
                        spe1_cover_cns = 'no'


                    if spe2_region in store_spe2_region_contain_cns_dic:
                        spe2_cover_cns = 'yes'
                    else:
                        spe2_cover_cns = 'no'

                    if spe2_region in store_acrregion_atlas_acr_celltype_dic:

                        rice_atlas_acr_dic = store_acrregion_atlas_acr_celltype_dic[spe2_region]

                        rice_atlas_acr_str = ';'.join(list(rice_atlas_acr_dic.keys()))

                        rice_atlas_celltye_list = []
                        for eachacr in rice_atlas_acr_dic:
                            rice_atlas_celltye_list.append(rice_atlas_acr_dic[eachacr])

                        rice_atlas_celltype_str = ';'.join(rice_atlas_celltye_list)

                    else:
                        rice_atlas_acr_str = 'none'
                        rice_atlas_celltype_str = 'none'

                    ##  store_eachpeak_celltype_str_num_dic[eachpeak] = {'celltypestr': store_celltype_str,
                     #                                            'celltypenum': str(celltypenum)}

                    ##check the cell type number
                    if rice_atlas_acr_str != 'none':


                        celltypenum_list = []
                        rice_atlas_acr_list = rice_atlas_acr_str.split(';')
                        for eachacr_atlas in rice_atlas_acr_list:
                            if eachacr_atlas in store_eachpeak_celltype_str_num_dic:
                                celltypenum = int(store_eachpeak_celltype_str_num_dic[eachacr_atlas]['celltypenum'])
                                celltypenum_list.append(celltypenum)

                        final_celltypenum = str(int(mean(celltypenum_list)))

                        ##check the number of organs
                        store_organ_dic = {}
                        for eachacr_atlas in rice_atlas_acr_list:
                            if eachacr_atlas in store_eachpeak_celltype_str_num_dic:
                                celltype_string = store_eachpeak_celltype_str_num_dic[eachacr_atlas]['celltypestr']
                                celltype_list = celltype_string.split(',')
                                for eachcelltype in celltype_list:
                                    mt = re.match('(.+)\.+',eachcelltype)
                                    organnm = mt.group(1)
                                    store_organ_dic[organnm] = 1

                        organnum = str(len(list(store_organ_dic.keys())))

                        ##updating 122423
                        ##if we find the cell type str all contains the leaf we will see there is no type cover
                        if len(list(store_organ_dic.keys())) == 1:

                            if list(store_organ_dic.keys())[0] == 'leaf':
                                organnum = 'none'
                                final_celltypenum = 'none'
                                rice_atlas_acr_str = 'none'


                    else:
                        organnum = 'none'
                        final_celltypenum = 'none'


                    final_line = eachline + '\t' + rice_atlas_acr_str + '\t' + rice_atlas_celltype_str + '\t' + spe1_cover_cns + '\t' + spe2_cover_cns + '\t' + organnum + '\t' + final_celltypenum
                    store_final_line_list.append(final_line)

                else:
                    final_line = eachline + '\t' + 'Atlas_ACR' + '\t' + 'Atlas_ACR_celltype' + '\t' + spe1_prefix + '_ACR_capture_CNS' + '\t' + spe2_prefix + '_region_capture_CNS' + '\t' + 'capture_atlas_organ_num' + '\t' + 'capture_altas_celltypenum'
                    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_to_' + spe2_prefix + '_summary_add_' + spe2_prefix + '_atlas_acr.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


def subfunction_add_check_num_celltypeACR (ipt_summary_fl,ipt_rice_atlas_fl,input_rice_atlas_acr_coverage_fl,
                                           s2_s2_organct_coverage_cutoff,spe1_prefix,spe2_prefix,opt_dir):

    ##updating 01112
    ##we will check the
    ##we will check the H3K27m3 for each of organ one by one

    ##updating 012424
    ##Here we will check directly use the ACRs from the atlas acr file
    store_all_rice_atlas_peak_dic = {}
    store_all_acr_str_dic = {}
    store_celltype_totalACR_dic = {}
    with open(ipt_rice_atlas_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            store_all_rice_atlas_peak_dic[acrnm] = 1

            acrstr = col[4]
            if acrstr == 'NA':
                new_acrstr = 'Broad'
            else:
                new_acrstr = acrstr

            store_all_acr_str_dic[acrnm] = new_acrstr

            new_acrstr_list = new_acrstr.split(',')
            for celltypestr in new_acrstr_list:
                if celltypestr in store_celltype_totalACR_dic:
                    store_celltype_totalACR_dic[celltypestr][acrnm] = 1
                else:
                    store_celltype_totalACR_dic[celltypestr] = {}
                    store_celltype_totalACR_dic[celltypestr][acrnm] = 1




    ##we still need the coverage plot information
    store_peak_organct_val_dic = {}
    store_total_organct_dic = {}
    with open(input_rice_atlas_acr_coverage_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            peaknm = col[1]
            organct = col[0]
            coverage_val = col[2]
            posneg_cate = col[3]
            if peaknm in store_peak_organct_val_dic:
                store_peak_organct_val_dic[peaknm][organct] = {'coverage_val': coverage_val,
                                                               'posneg_cate': posneg_cate}
            else:
                store_peak_organct_val_dic[peaknm] = {}
                store_peak_organct_val_dic[peaknm][organct] = {'coverage_val': coverage_val,
                                                               'posneg_cate': posneg_cate}

            store_total_organct_dic[organct] = 1

    store_final_line_list = []
    store_eachpeak_celltype_str_num_dic = {}
    #store_celltype_totalACR_dic = {}
    for eachpeak in store_peak_organct_val_dic:

        ##we want to allow the peaks are in the final peaks we called
        if eachpeak in store_all_rice_atlas_peak_dic:

            store_celltype_list = []
            store_cover_val_list = []

            for eachorganct in store_peak_organct_val_dic[eachpeak]:
                coverage_val = store_peak_organct_val_dic[eachpeak][eachorganct]['coverage_val']
                posneg_cate = store_peak_organct_val_dic[eachpeak][eachorganct]['posneg_cate']

                if float(coverage_val) >= float(s2_s2_organct_coverage_cutoff):

                    ##updating 050423
                    if posneg_cate == 'Positive':
                        store_celltype_list.append(eachorganct)
                        store_cover_val_list.append(coverage_val)

            store_celltype_str = ','.join(store_celltype_list)
            store_cover_val_str = ','.join(store_cover_val_list)

            if store_celltype_str != '':
                celltypenum = len(store_celltype_list)
                final_line = 'CutoffCover:' + s2_s2_organct_coverage_cutoff + '\t' + eachpeak + '\t' + store_celltype_str + '\t' + store_cover_val_str + '\t' + str(
                    celltypenum)
                store_final_line_list.append(final_line)
                store_eachpeak_celltype_str_num_dic[eachpeak] = {'celltypestr': store_celltype_str,
                                                                 'celltypenum': str(celltypenum)}

                #celltype_list = store_celltype_str.split(',')

                #for eachcelltype in celltype_list:
                #    if eachcelltype in store_celltype_totalACR_dic:
                #        store_celltype_totalACR_dic[eachcelltype][eachpeak] = 1
                #    else:
                #        store_celltype_totalACR_dic[eachcelltype] = {}
                #        store_celltype_totalACR_dic[eachcelltype][eachpeak] = 1













    store_celltype_ACR_dic = {}
    store_celltype_ACR_not_coverCNS_dic = {}
    store_final_line_rmLeafCT_list = []
    count = 0
    ##the example is the maize to rice
    with open (ipt_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split()
                atlasACR = col[19]
                riceACR = col[5]
                maizeACRcaptureCNS = col[21]

                ##updating 020224
                riceregioncaptureCNS = col[22]

                maizeID = col[0]


                if atlasACR in store_eachpeak_celltype_str_num_dic:
                    celltypestr = store_eachpeak_celltype_str_num_dic[atlasACR]['celltypestr']
                    celltypelist = celltypestr.split(',')
                    otherorganct_dic = {}
                    for eachcelltypestr in celltypelist:

                        mt = re.match('(.+)\..+',eachcelltypestr)
                        organ = mt.group(1)
                        if organ != 'leaf':
                            otherorganct_dic[eachcelltypestr] = 1


                    new_updated_celltypenum = str(len(list(otherorganct_dic.keys())))

                else:
                    new_updated_celltypenum = 'none'

                ##we will select the variable ACR other than all the ACR within the syntenic regions
                if atlasACR != 'none' and riceACR == 'none' and maizeID != 'none':

                    if maizeACRcaptureCNS == 'yes' :

                        if riceregioncaptureCNS == 'yes':

                            atlasacr_list = atlasACR.split(';')

                            for eachatlasACR in atlasacr_list:

                                ##Here we will check the celltype str
                                celltypestr = store_all_acr_str_dic[eachatlasACR]

                                celltype_list = celltypestr.split(',')

                                for eachcelltype in celltype_list:

                                    if eachcelltype in store_celltype_ACR_dic:
                                        store_celltype_ACR_dic[eachcelltype][eachatlasACR] = 1
                                    else:
                                        store_celltype_ACR_dic[eachcelltype] = {}
                                        store_celltype_ACR_dic[eachcelltype][eachatlasACR] = 1

                    else:

                        atlasacr_list = atlasACR.split(';')

                        for eachatlasACR in atlasacr_list:

                            ##Here we will check the celltype str
                            celltypestr = store_all_acr_str_dic[eachatlasACR]

                            celltype_list = celltypestr.split(',')

                            for eachcelltype in celltype_list:

                                if eachcelltype in store_celltype_ACR_not_coverCNS_dic:
                                    store_celltype_ACR_not_coverCNS_dic[eachcelltype][eachatlasACR] = 1
                                else:
                                    store_celltype_ACR_not_coverCNS_dic[eachcelltype] = {}
                                    store_celltype_ACR_not_coverCNS_dic[eachcelltype][eachatlasACR] = 1


                final_line = eachline + '\t' + new_updated_celltypenum
                store_final_line_rmLeafCT_list.append(final_line)

            else:
                final_line = eachline + '\t' + 'capture_altas_celltypenum_noleaf'
                store_final_line_rmLeafCT_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_to_' + spe2_prefix + '_summary_add_rice_atlas_acr_rmleaf.txt','w+') as opt:
        for eachline in store_final_line_rmLeafCT_list:
            opt.write(eachline + '\n')




    store_final_line_list = []
    for eachcelltype in store_celltype_ACR_dic:

        celltype_CNS_atlasACR_count = len(list(store_celltype_ACR_dic[eachcelltype].keys()))

        celltype_total_atlasACR_count = len(list(store_celltype_totalACR_dic[eachcelltype].keys()))

        final_line = 'OverlapCNS' + '\t' + eachcelltype + '\t' + str(celltype_CNS_atlasACR_count) + '\t' + str(celltype_total_atlasACR_count) + '\t' + \
                     str(celltype_CNS_atlasACR_count/celltype_total_atlasACR_count)

        store_final_line_list.append(final_line)



    for eachcelltype in store_celltype_ACR_not_coverCNS_dic:

        celltype_CNS_atlasACR_count = len(list(store_celltype_ACR_not_coverCNS_dic[eachcelltype].keys()))

        celltype_total_atlasACR_count = len(list(store_celltype_totalACR_dic[eachcelltype].keys()))

        final_line = 'NotOverlapCNS' + '\t' + eachcelltype + '\t' + str(celltype_CNS_atlasACR_count) + '\t' + str(celltype_total_atlasACR_count) + '\t' + \
                     str(celltype_CNS_atlasACR_count/celltype_total_atlasACR_count)

        store_final_line_list.append(final_line)


    with open (opt_dir + '/opt_celltype_ACRcount_CNSaltasACR_' + spe1_prefix + '_' + spe2_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



    store_all_celltype_dic = {}
    for eachcelltype in store_celltype_ACR_dic:
        store_all_celltype_dic[eachcelltype] = 1
    for eachcelltype in store_celltype_ACR_not_coverCNS_dic:
        store_all_celltype_dic[eachcelltype] = 1

    store_check_relative_prop_line_list = []
    first_line = 'Celltype' + '\t'  + 'FinalCelltype' + '\t' + 'Organ' +  '\t' + 'OverlapCNS_ACRcount' + '\t' + 'NotOverlapCNS_ACRcount' + '\t' + 'OverlapCNS_ACRprop'
    store_check_relative_prop_line_list.append(first_line)
    for eachcelltype in store_all_celltype_dic:

        if eachcelltype in store_celltype_ACR_dic:
            celltype_CNS_atlasACR_overlapCNS_count = len(list(store_celltype_ACR_dic[eachcelltype].keys()))
        else:
            celltype_CNS_atlasACR_overlapCNS_count = 0

        if eachcelltype in store_celltype_ACR_not_coverCNS_dic:
            celltype_CNS_atlasACR_notoverlapCNS_count = len(list(store_celltype_ACR_not_coverCNS_dic[eachcelltype].keys()))
        else:
            celltype_CNS_atlasACR_notoverlapCNS_count = 0

        #count_list = [celltype_CNS_atlasACR_overlapCNS_count,celltype_CNS_atlasACR_notoverlapCNS_count]
        #max_count = max(count_list)

        #abs_diff = abs(celltype_CNS_atlasACR_overlapCNS_count - celltype_CNS_atlasACR_notoverlapCNS_count)

        #rel_difference_prop = abs_diff/max_count

        total_count = celltype_CNS_atlasACR_overlapCNS_count + celltype_CNS_atlasACR_notoverlapCNS_count
        overlapCNS_prop = celltype_CNS_atlasACR_overlapCNS_count/total_count

        ##make more categories
        if re.match('(.+)\.(.+)',eachcelltype):
            mt = re.match('(.+)\.(.+)',eachcelltype)
            organ = mt.group(1)
            celltype = mt.group(2)

            ##for the leaf related cell types
            final_celltype = 'notInclude'
            if organ == 'bud' or organ == 'Eseedling' or organ == 'seedling':

                if 'unknown' not in celltype and 'Bulliform' != celltype and \
                        'LargeParenchyma' != 'celltype' and 'VascularRelateCell' != celltype and 'SAM' != celltype and 'VascularSclerenchyma' != celltype:

                    if 'DevelopingFiber' == celltype or 'Epidermis' == celltype or 'GuardCell' == celltype or \
                         'StomatalPrecursor' == celltype:

                        final_celltype = 'EpidermalRelatedCell'

                    else:
                        if 'Mesophyll' == celltype:
                            final_celltype = 'Mesophyll'
                        else:
                            if 'BundleSheath' == celltype:
                                final_celltype = 'BundleSheath'
                            else:

                                if 'Parenchyma' in celltype:

                                    final_celltype = 'VasParenchyma'

                                else:

                                    final_celltype = 'CCandVasElement'

            if organ == 'semroot' or organ == 'crownroot':

                if 'unknown' not in celltype:

                    if celltype == 'QC':
                        final_celltype = 'QC'
                    else:
                        if celltype == 'Columella':
                            final_celltype = 'Columella'
                        else:
                            if celltype == 'Atrichoblast':
                                final_celltype = 'Atrichoblast'
                            else:
                                if celltype == 'Exodermis':
                                    final_celltype = 'Exodermis'
                                else:
                                    if celltype == 'Endodermis' or celltype == 'Cortex' or celltype == 'CortexMeristem':
                                        final_celltype = 'GroundTissueCell'
                                    else:
                                        final_celltype = 'SteleCell'

            if organ == 'Lseed' or organ == 'Eseed':

                if 'unknown' not in celltype:

                    if celltype == 'CSE' or celltype == 'DSE' or celltype == 'AL' or celltype == 'DSELSE' or celltype == 'LSE':

                        final_celltype = 'EndospermCell'

                    else:
                        if celltype == 'SA' or celltype == 'Epiblast' or celltype == 'NucellarEpidermis' or celltype == 'Coleoptile' or \
                            celltype == 'SVB' or celltype == 'Scutellum' or celltype == 'RA' or celltype == 'Coleorhiza':

                            final_celltype = 'Embryo'

            if organ == 'panicle':

                if 'unknown' not in celltype:

                    if celltype == 'SM' or celltype == 'SBM' or celltype == 'PBM' or celltype == 'FM':
                        final_celltype = 'MeristemCell'
                    else:
                        if celltype == 'EpidermisMeri':
                            final_celltype = 'EpidermisMeri'
                        else:
                            if celltype != 'Rachis' and celltype != 'IM' and celltype != 'Bract':
                                final_celltype = 'FloretCell'

        else:

            final_celltype = 'Broad'
            organ = 'All'

        final_line = eachcelltype + '\t' + final_celltype + '\t' + organ + '\t' + \
                     str(celltype_CNS_atlasACR_overlapCNS_count) + '\t' + str(celltype_CNS_atlasACR_notoverlapCNS_count) + '\t' + str(overlapCNS_prop)
        store_check_relative_prop_line_list.append(final_line)



    with open (opt_dir + '/opt_celltype_ACRcount_CNSaltasACR_' + spe1_prefix + '_' + spe2_prefix + '_relativeDiffProp.txt','w+') as opt:
        for eachline in store_check_relative_prop_line_list:
            opt.write(eachline + '\n')




















