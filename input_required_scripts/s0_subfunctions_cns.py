#!/usr/bin/env python

##this script is to check the basics of cns overlapping each of species ACRs

import re
import glob
import sys
import subprocess
import os


def subfunction_check_cns_overlap_ACR (ipt_cell_type_spe_fl,ipt_spe_cns_fl, ipt_gene_cate_fl,opt_dir,
                                       spe_prefix):


    ##Here it would be extracting the cns file first
    number_col = 0
    with open(ipt_spe_cns_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split('\t')
                number_col = len(col)

    if number_col > 4:

        store_final_line_list = []
        with open(ipt_spe_cns_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                if not eachline.startswith('#'):
                    col = eachline.strip().split()
                    cns_line = col[0] + '\t' + col[3] + '\t' + col[4]
                    store_final_line_list.append(cns_line)

    else:
        store_final_line_list = []
        with open(ipt_spe_cns_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                final_line = col[0] + '\t' + col[1] + '\t' + col[2]
                store_final_line_list.append(final_line)

    with open(opt_dir + '/temp_cns.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_cns.txt > ' + \
          opt_dir + '/temp_cns_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    cmd = 'cut -f 1-3 ' + ipt_cell_type_spe_fl + ' > ' + opt_dir + '/temp_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_acr.txt > ' + opt_dir + '/temp_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##now we will do the intersecting to find the CNS
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_acr_sorted.txt' + ' -b ' + opt_dir + '/temp_cns_sorted.txt' + ' > ' + \
          opt_dir + '/temp_' + spe_prefix + '_celltype_cns_intersect.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)


    store_ACR_cns_num_dic = {}
    with open (opt_dir + '/temp_' + spe_prefix + '_celltype_cns_intersect.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            cnsloc = col[3] + '_' + col[4] + '_' + col[5]

            if acrloc in store_ACR_cns_num_dic:
                store_ACR_cns_num_dic[acrloc][cnsloc] = 1
            else:
                store_ACR_cns_num_dic[acrloc] = {}
                store_ACR_cns_num_dic[acrloc][cnsloc] = 1

    store_ACR_cate_dic = {}
    with open (ipt_gene_cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0]
            genecate = col[1]
            store_ACR_cate_dic[acrloc] = genecate

    store_final_line_list = []
    with open(ipt_cell_type_spe_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            mt = re.match('.+;(.+)',col[3])
            celltype = mt.group(1)

            if acrloc in store_ACR_cns_num_dic:
                cnsloc_num = len(list(store_ACR_cns_num_dic[acrloc].keys()))
            else:
                cnsloc_num = '0'

            if acrloc in store_ACR_cate_dic:
                genecate = store_ACR_cate_dic[acrloc]
            else:
                genecate = 'none'

            final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + celltype + '\t' + str(cnsloc_num) + '\t' + genecate
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe_prefix + '_celltype_cns_num.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##updating 122323
    ##we will use the cns as the subject
    store_cns_overlapACR_dic = {}
    with open(opt_dir + '/temp_' + spe_prefix + '_celltype_cns_intersect.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            cnsloc = col[3] + '_' + col[4] + '_' + col[5]

            if cnsloc in store_cns_overlapACR_dic:
                store_cns_overlapACR_dic[cnsloc][acrloc] = 1
            else:
                store_cns_overlapACR_dic[cnsloc] = {}
                store_cns_overlapACR_dic[cnsloc][acrloc] = 1

    store_acr_celltype_dic = {}
    with open(ipt_cell_type_spe_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            mt = re.match('.+;(.+)', col[3])
            celltype = mt.group(1)
            store_acr_celltype_dic[acrloc] = celltype


    store_final_line_list = []
    with open (opt_dir + '/temp_cns_sorted.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            cnsloc = col[0] + '_' + col[1] + '_' + col[2]
            if cnsloc in store_cns_overlapACR_dic:
                acrloc_dic = store_cns_overlapACR_dic[cnsloc]

                acrloc_list = list(acrloc_dic.keys())
                celltype_list = []
                for eachacrloc in acrloc_dic:

                    if eachacrloc in store_acr_celltype_dic:
                        celltype = store_acr_celltype_dic[eachacrloc]
                    else:
                        celltype = 'none'

                    celltype_list.append(celltype)

                acrloc_str = ','.join(acrloc_list)
                celltype_str = ','.join(celltype_list)

                ##generate a final one
                if len(celltype_list) != 1:

                    broad_acc_count = 0
                    for eachcellytpe in celltype_list:
                        if 'broadly_accessible' == eachcellytpe:
                            broad_acc_count += 1

                    if broad_acc_count != len(celltype_list):
                        final_celltype_cate = 'CT'
                    else:
                        final_celltype_cate = 'Broad'

                else:
                    if 'broadly_accessible' == celltype_list[0]:
                        final_celltype_cate = 'Broad'
                    else:
                        final_celltype_cate = 'CT'


            else:
                acrloc_str = 'none'
                celltype_str = 'none'
                final_celltype_cate = 'none'

            final_line = eachline + '\t' + acrloc_str + '\t' + celltype_str + '\t' + final_celltype_cate
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe_prefix + '_cns_subject_acr_celltype.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')





def subfunction_intersect_with_atlas (ipt_cell_type_rice_fl,ipt_atlas_acr_fl,ipt_spe_cns_fl,opt_dir):

    number_col = 0
    with open(ipt_spe_cns_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split('\t')
                number_col = len(col)

    if number_col > 4:

        store_final_line_list = []
        with open(ipt_spe_cns_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                if not eachline.startswith('#'):
                    col = eachline.strip().split()
                    cns_line = col[0] + '\t' + col[3] + '\t' + col[4]
                    store_final_line_list.append(cns_line)

    else:
        store_final_line_list = []
        with open(ipt_spe_cns_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                final_line = col[0] + '\t' + col[1] + '\t' + col[2]
                store_final_line_list.append(final_line)

    with open(opt_dir + '/temp_cns.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_cns.txt > ' + \
          opt_dir + '/temp_cns_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##intersect the atlas to the cell type spe fl
    cmd = 'cut -f 1-3 ' + ipt_cell_type_rice_fl + ' > ' + opt_dir + '/temp_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_acr.txt > ' + opt_dir + '/temp_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + ipt_atlas_acr_fl + ' -b ' + opt_dir + '/temp_acr_sorted.txt'  + ' > ' + opt_dir + '/temp_intersect_rice_altas_leafACR.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_rice_atlas_acr_leafACR_dic = {}
    with open (opt_dir + '/temp_intersect_rice_altas_leafACR.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            leafacrnm = col[5] + '_' + col[6] + '_' + col[7]
            if acrnm in store_rice_atlas_acr_leafACR_dic:
                store_rice_atlas_acr_leafACR_dic[acrnm][leafacrnm] = 1
            else:
                store_rice_atlas_acr_leafACR_dic[acrnm] = {}
                store_rice_atlas_acr_leafACR_dic[acrnm][leafacrnm] = 1

    ##intersect the atlas ACR to spe_cns data
    cmd = 'bedtools intersect -wa -wb -a ' + ipt_atlas_acr_fl + ' -b ' + opt_dir + '/temp_cns_sorted.txt > ' + \
          opt_dir + '/temp_intersect_rice_atlas_cns.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_rice_atlas_cnsloc_dic = {}
    with open (opt_dir + '/temp_intersect_rice_atlas_cns.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            cnsloc = col[5] + '_' + col[6] + '_' + col[7]

            if acrnm in store_rice_atlas_cnsloc_dic:
                store_rice_atlas_cnsloc_dic[acrnm][cnsloc] = 1
            else:
                store_rice_atlas_cnsloc_dic[acrnm] = {}
                store_rice_atlas_cnsloc_dic[acrnm][cnsloc] = 1

    store_final_line_list = []
    with open (ipt_atlas_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            acrnm = col[0] + '_' + col[1] + '_' + col[2]

            if acrnm in store_rice_atlas_acr_leafACR_dic:
                leafACR_str = ','.join(list(store_rice_atlas_acr_leafACR_dic[acrnm].keys()))
            else:
                leafACR_str = 'none'

            if acrnm in store_rice_atlas_cnsloc_dic:
                cns_num = str(len(list(store_rice_atlas_cnsloc_dic[acrnm].keys())))
            else:
                cns_num = '0'

            final_line = eachline + '\t' + leafACR_str + '\t' + cns_num
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_rice_atlas_acr_addPabloLeaf_addcns.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')






















