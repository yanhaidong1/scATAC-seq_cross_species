#!/usr/bin/env python

##updating 120223 we will add the log10tn5 and library
##Here we will check the TF accessiblity

import re
import glob
import sys
import subprocess
import os
import random
import math

##we will first check number of cells per cell types
##decide the number of cells we will compare for the four cell types
##just focus on four cell types



def subfunction_flt_meta_cells_keep_same_cellnum_celltypes (ipt_meta_dir,s0_s6_sub_cluster_colnum_str,opt_dir):


    store_spe_clusternum_dic = {}
    for eachpair in s0_s6_sub_cluster_colnum_str.split(','):
        mt = re.match('(.+):(.+)',eachpair)
        store_spe_clusternum_dic[mt.group(1)] = mt.group(2)

    all_meta_fl_list = glob.glob(ipt_meta_dir + '/*')

    store_spe_cellcluster_num_dic = {}
    store_spe_cellcluster_celllist_dic = {}
    for eachmetafl in all_meta_fl_list:

        mt = re.match('.+/(.+)',eachmetafl)
        flnm = mt.group(1)

        mt = re.match('(.+)_meta\.txt',flnm)
        spenm = mt.group(1)

        store_cellcluster_celllist_dic = {}
        count = 0
        with open (eachmetafl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                cellnm = col[0]
                count += 1
                if count != 1:
                    cellcluster = col[int(store_spe_clusternum_dic[spenm]) - 1]
                    #print(cellcluster)

                    if cellcluster == 'bundle_sheath' or \
                            cellcluster == 'epidermis' or \
                            cellcluster == 'mesophyll' or \
                            cellcluster == 'companion_cells_sieve_elements' or \
                            cellcluster == 'companion_cell' or \
                            cellcluster == 'bundle_sheath_ncell_7219' or \
                            cellcluster == 'companion_cells_sieve_elements_ncell_1193' or \
                            cellcluster == 'epidermis_ncell_5440' or \
                            cellcluster == 'mesophyll_ncell_5257':

                        spe_cellclustnm = spenm + '_' + cellcluster
                        if spe_cellclustnm in store_spe_cellcluster_num_dic:
                            store_spe_cellcluster_num_dic[spe_cellclustnm] += 1
                        else:
                            store_spe_cellcluster_num_dic[spe_cellclustnm] = 1

                        if spe_cellclustnm in store_cellcluster_celllist_dic:
                            store_cellcluster_celllist_dic[spe_cellclustnm].append(cellnm)
                        else:
                            store_cellcluster_celllist_dic[spe_cellclustnm] = []
                            store_cellcluster_celllist_dic[spe_cellclustnm].append(cellnm)

        store_spe_cellcluster_celllist_dic[spenm] = store_cellcluster_celllist_dic

    ##find the least number
    least_num_key = min(store_spe_cellcluster_num_dic, key=store_spe_cellcluster_num_dic.get)
    least_num = store_spe_cellcluster_num_dic[least_num_key]

    ##randomly pick cells per cell type
    store_spe_select_cellcluster_celllist_dic = {}
    for eachspe in store_spe_cellcluster_celllist_dic:
        store_cellcluster_celllist_dic = store_spe_cellcluster_celllist_dic[eachspe]

        store_cellcluster_celllist = {}
        for eachcellcluster in store_cellcluster_celllist_dic:
            celllist = store_cellcluster_celllist_dic[eachcellcluster]

            rd_celllist = random.sample(celllist, least_num)

            store_cellcluster_celllist[eachcellcluster] = rd_celllist

        store_spe_select_cellcluster_celllist_dic[eachspe] = store_cellcluster_celllist

    #print(store_spe_select_cellcluster_celllist_dic['maize']['maize_protoderm'])

    ##generate a new meta file with least num
    for eachmetafl in all_meta_fl_list:

        mt = re.match('.+/(.+)',eachmetafl)
        flnm = mt.group(1)

        mt = re.match('(.+)_meta\.txt',flnm)
        spenm = mt.group(1)

        print(spenm)

        store_final_line_list = []
        count = 0
        with open (eachmetafl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                cellnm = col[0]
                count += 1
                if count != 1:
                    cellcluster = col[int(store_spe_clusternum_dic[spenm]) - 1]

                    if cellcluster == 'bundle_sheath' or \
                            cellcluster == 'epidermis' or \
                            cellcluster == 'mesophyll' or \
                            cellcluster == 'companion_cells_sieve_elements' or \
                            cellcluster == 'companion_cell' or \
                            cellcluster == 'bundle_sheath_ncell_7219' or \
                            cellcluster == 'companion_cells_sieve_elements_ncell_1193' or \
                            cellcluster == 'epidermis_ncell_5440' or \
                            cellcluster == 'mesophyll_ncell_5257':

                        print(cellcluster)
                        spe_clusternm = spenm + '_' + cellcluster
                        filter_cell_list = store_spe_select_cellcluster_celllist_dic[spenm][spe_clusternm]

                        if cellnm in filter_cell_list:
                            store_final_line_list.append(eachline)

                else:
                    store_final_line_list.append(eachline)

        with open (opt_dir + '/' + spenm + '_flt_meta.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')



def subfunction_build_peak_cell_in_meta (input_meta_fl,input_peak_cell_sparse_fl,input_output_dir,spe_prefix):

    store_cell_peak_num_dic = {}
    with open (input_peak_cell_sparse_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            peaknm = col[0]
            cell_nm = col[1]

            if peaknm.startswith('chr'):
                mt = re.match('chr(.+)', peaknm)
                peak_nm = mt.group(1)
            else:
                peak_nm = peaknm

            if cell_nm in store_cell_peak_num_dic:

                if peak_nm in store_cell_peak_num_dic[cell_nm]:
                    store_cell_peak_num_dic[cell_nm][peak_nm] += 1
                else:
                    store_cell_peak_num_dic[cell_nm][peak_nm] = 1

            else:
                store_cell_peak_num_dic[cell_nm] = {}
                store_cell_peak_num_dic[cell_nm][peak_nm] = 1


    store_final_line_list = []
    count = 0
    with open (input_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            count += 1
            if count != 1:
                cellnm = col[0]

                peak_dic = store_cell_peak_num_dic[cellnm]
                ACR_num = len(list(peak_dic.keys()))

                final_line = eachline + '\t' + str(ACR_num)
                store_final_line_list.append(final_line)

            else:
                final_line = eachline + '\t' + 'ACRnum'
                store_final_line_list.append(final_line)

    with open (input_output_dir + '/opt_' + spe_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


def subfunction_add_log10tn5_library (ipt_meta_fl,s0_s6_sub_target_library_colnum,s0_s6_sub_target_tn5_colnum,spe_prefix,
                                      opt_dir):

    store_final_line_list = []
    count = 0
    with open (ipt_meta_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                library_col = col[int(s0_s6_sub_target_library_colnum) - 1]
                tn5_colnm = col[int(s0_s6_sub_target_tn5_colnum) - 1]

                final_line = eachline + '\t' + str(math.log10(int(tn5_colnm))) + '\t' + library_col
                store_final_line_list.append(final_line)
            else:
                final_line = eachline + '\t' + 'log10tn5' + '\t' + 'library'
                store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



def subfunction_intersect_acr_motif (ipt_spe_ACR_fl,ipt_motif_fl,opt_dir,spe_prefix):

    cmd = 'bedtools intersect -wa -wb -a ' + ipt_spe_ACR_fl + ' -b ' + ipt_motif_fl + ' > ' + \
          opt_dir + '/opt_' + spe_prefix + '_ACR_motif.bed'
    subprocess.call(cmd,shell=True)

    store_total_motif_nm_dic = {}
    store_final_line_list = []
    with open(opt_dir + '/opt_' + spe_prefix + '_ACR_motif.bed', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            if col[8].startswith('MA'):
                motif_nm = col[8]
            else:
                motif_nm = col[9]

            mt = re.match('(.+)_.+', motif_nm)
            concise_motif = mt.group(1)

            final_line = col[0]
            for i in range(1, 8):
                final_line = final_line + '\t' + col[i]

            final_line = final_line + '\t' + concise_motif + '\t' + col[9] + '\t' + col[10] + '\t' + col[11]
            store_final_line_list.append(final_line)

            store_total_motif_nm_dic[concise_motif] = 1

    with open(opt_dir + '/opt_' + spe_prefix + '_ACR_motif_modinm.bed', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    with open(opt_dir + '/opt_' + spe_prefix + '_total_motif_nm.txt', 'w+') as opt:
        for eachline in store_total_motif_nm_dic:
            opt.write(eachline + '\n')


def subfunction_prepare_regr_enrich_motif (input_required_script_dir,ipt_acr_motif_modinm_fl,
                                           ipt_peak_sparse_fl,s0_s6_sub_pvalcutoff,spe_prefix,
                                           opt_dir
                                           ):

    print('open step01_s6_prepare_regr_enrich_fimoMotif')

    s6_script_opt2 = input_required_script_dir + \
                     '/s0_subfunctions_motif_model_enrich/motif_enrichment/prepare_regression_data.R'

    #motif_acr_fl = s2_intersect_acr_dir + '/opt_ACR_motif_modinm.bed'
    motif_acr_fl = ipt_acr_motif_modinm_fl

    cmd = 'Rscript ' + s6_script_opt2 + \
          ' ' + motif_acr_fl + \
          ' ' + s0_s6_sub_pvalcutoff + \
          ' ' + ipt_peak_sparse_fl + \
          ' ' + spe_prefix + \
          ' ' + opt_dir
    print(cmd)
    subprocess.call(cmd, shell=True)




def subfunction_conduct_enrichment (input_required_script_dir,ipt_meta_add_ACRnum_fl,
                                    ipt_prepare_regr_enrich_dir_dir,opt_dir,s0_s6_sub_target_clusternm,
                                    input_core_num
                                    ):

    s7_script_opt2 = input_required_script_dir + \
                     '/s0_subfunctions_motif_model_enrich/motif_enrichment/conduct_enrichment.R'

    #target_meta_fl = s5_prepare_regr_enrich_dir + '/opt_' + step01_s5_prefix_meta + '.txt'

    cmd = 'Rscript ' + s7_script_opt2 + \
          ' ' + ipt_prepare_regr_enrich_dir_dir + \
          ' ' + ipt_meta_add_ACRnum_fl + \
          ' ' + opt_dir + \
          ' ' + s0_s6_sub_target_clusternm + \
          ' ' + input_core_num
    print(cmd)
    subprocess.call(cmd, shell=True)
