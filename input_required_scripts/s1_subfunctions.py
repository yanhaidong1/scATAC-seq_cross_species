#!/usr/bin/env python



import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats

def subfunction_iden_correspond (spe1_spe2_blast_fl,spe1_spe2_ref_blast_fl,spe1_H3K27me3_fl,spe2_H3K27me3_fl,
                                 spe1_acr_fl,spe1_spe2_syntenic_acr_fl,opt_dir,spe1_prefix,spe2_prefix):

    cmd = 'cut -f 1-3 ' + spe2_H3K27me3_fl + ' > ' + opt_dir + '/temp_' + spe2_prefix + '_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe2_prefix + '_acr.txt > ' + opt_dir + '/temp_' + spe2_prefix + '_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_spe1_ID_spe1_acr_dic = {}
    store_spe1_acr_spe1_ID_dic = {}
    with open (spe1_spe2_ref_blast_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            acrID = col[3].split(';')[0]
            store_spe1_ID_spe1_acr_dic[acrID] = acrnm
            store_spe1_acr_spe1_ID_dic[acrnm] = acrID

    store_spe1_acrID_spe2_region_dic = {}
    store_final_line_list = []
    with open (spe1_spe2_blast_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            spe1_acrID = col[3].split(';')[1]
            spe2_region = col[0] + '_' + col[1] + '_' + col[2]
            if spe1_acrID in store_spe1_acrID_spe2_region_dic:
                store_spe1_acrID_spe2_region_dic[spe1_acrID][spe2_region] = 1
            else:
                store_spe1_acrID_spe2_region_dic[spe1_acrID] = {}
                store_spe1_acrID_spe2_region_dic[spe1_acrID][spe2_region] = 1

            final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[0] + '_' + col[1] + '_' + col[2]
            store_final_line_list.append(final_line)

    with open (opt_dir + '/temp_' + spe2_prefix + '_region_loc.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe2_prefix + '_region_loc.txt > ' + opt_dir + '/temp_' + spe2_prefix + '_region_loc_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + spe2_prefix + '_region_loc_sorted.txt' + ' -b ' + \
          opt_dir + '/temp_' + spe2_prefix + '_acr_sorted.txt > ' + opt_dir + '/temp_intersect_' + spe2_prefix + '_region_ACR.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_spe2_region_acr_dic = {}
    with open (opt_dir + '/temp_intersect_' + spe2_prefix + '_region_ACR.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            region_loc = col[3]
            acr = col[4] + '_' + col[5] + '_' + col[6]
            if region_loc in store_spe2_region_acr_dic:
                store_spe2_region_acr_dic[region_loc][acr] = 1
            else:
                store_spe2_region_acr_dic[region_loc] = {}
                store_spe2_region_acr_dic[region_loc][acr] = 1

    ##for the H3K27me3 file
    store_target_ACR_cate_spe1_dic = {}
    count = 0
    with open(spe1_H3K27me3_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            broad_restrict_cate = col[3]
            final_cate = col[-1]
            count += 1
            if count != 1:
                store_target_ACR_cate_spe1_dic[acrnm] = {'BroOrRes': broad_restrict_cate, 'InorNotH3K27m3': final_cate}

    store_target_ACR_cate_spe2_dic = {}
    store_final_line_list = []
    count = 0
    with open(spe2_H3K27me3_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            broad_restrict_cate = col[3]
            final_cate = col[-1]
            count += 1

            if spe2_prefix != 'rice':

                if count != 1:
                    store_target_ACR_cate_spe2_dic[acrnm] = {'BroOrRes': broad_restrict_cate, 'InorNotH3K27m3': final_cate}
                    acrline = col[0] + '\t' + col[1] + '\t' + col[2]
                    store_final_line_list.append(acrline)

            else:
                store_target_ACR_cate_spe2_dic[acrnm] = {'BroOrRes': broad_restrict_cate, 'InorNotH3K27m3': final_cate}
                acrline = col[0] + '\t' + col[1] + '\t' + col[2]
                store_final_line_list.append(acrline)


    ##updating 101723
    ##we will store the cell type specific ACRs
    store_spe1_ACR_celltypecate_dic = {}
    with open (spe1_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split ()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            mt = re.match('(.+);(.+)',col[3])
            acrID = mt.group(1)
            celltypecate = mt.group(2)
            store_spe1_ACR_celltypecate_dic[acrnm] = celltypecate

    store_spe1_ACR_syntenic_region_dic =  {}
    with open (spe1_spe2_syntenic_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            syntenicID = col[0]
            ACR_str = col[1]
            ACR_list = ACR_str.split(',')
            for eachACR in ACR_list:
                if eachACR in store_spe1_ACR_syntenic_region_dic:
                    store_spe1_ACR_syntenic_region_dic[eachACR].append(syntenicID)
                else:
                    store_spe1_ACR_syntenic_region_dic[eachACR] = []
                    store_spe1_ACR_syntenic_region_dic[eachACR].append(syntenicID)


    ##we will do the following stuffs
    ##1) check how many spe1 ACR could have blast to the syntenic regions
    ##2) check how many spe1 ACR could overlap with the ACRs in spe2
    ##3) check how many syntenic_ACR_spe2_ACR having H3K27me3 and how many of them are broad or restricted
    ##4) check how many syntenic_ACR_spe2_region having H3K27me3 broad and how many of them are broad or restricted

    ##so the final file we will built is the
    ##spe1_ACR_loc   spe2_region_loc   spe2_acr  spe2_H3K27me3_overlapOrNot spe2_H3K27me3_overlapCate

    store_final_line_list = []
    first_line = spe1_prefix + '_ACRID' + '\t' + spe1_prefix + '_ACR' + '\t' + spe1_prefix + '_InorNotH3K27m3' + '\t' + spe1_prefix + '_BroOrRes' + '\t' + \
                 spe2_prefix + '_region' + '\t' + spe2_prefix + '_ACR' + '\t' + spe2_prefix + '_InorNotH3K27m3' + '\t' + spe2_prefix + '_BroOrRes'
    store_final_line_list.append(first_line)
    for eachspe1ACRID in store_spe1_acrID_spe2_region_dic:

        eachspe1ACR = store_spe1_ID_spe1_acr_dic[eachspe1ACRID]

        if eachspe1ACR in store_target_ACR_cate_spe1_dic:

            InorNotH3K27m3_spe1 = store_target_ACR_cate_spe1_dic[eachspe1ACR]['InorNotH3K27m3']
            BroOrRes_spe1 = store_target_ACR_cate_spe1_dic[eachspe1ACR]['BroOrRes']

            spe2_region_dic = store_spe1_acrID_spe2_region_dic[eachspe1ACRID]

            for eachspe2region in spe2_region_dic:

                if eachspe2region in store_spe2_region_acr_dic:
                    spe2_acr_dic = store_spe2_region_acr_dic[eachspe2region]

                    for eachspe2ACR in spe2_acr_dic:
                        InorNotH3K27m3_spe2 = store_target_ACR_cate_spe2_dic[eachspe2ACR]['InorNotH3K27m3']
                        BroOrRes_spe2 = store_target_ACR_cate_spe2_dic[eachspe2ACR]['BroOrRes']

                        final_line = eachspe1ACRID + '\t' + eachspe1ACR + '\t' + InorNotH3K27m3_spe1 + '\t' + BroOrRes_spe1 + '\t' + \
                                     eachspe2region + '\t' + eachspe2ACR + '\t' + InorNotH3K27m3_spe2 + '\t' + BroOrRes_spe2
                        store_final_line_list.append(final_line)

                else:
                    InorNotH3K27m3_spe2 = 'none'
                    BroOrRes_spe2 = 'none'
                    spe2_acr = 'none'

                    final_line = eachspe1ACRID + '\t' + eachspe1ACR + '\t' + InorNotH3K27m3_spe1 + '\t' + BroOrRes_spe1 + '\t' + \
                                 eachspe2region + '\t' + spe2_acr + '\t' + InorNotH3K27m3_spe2 + '\t' + BroOrRes_spe2
                    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##updating 101723
    store_final_line_list = []
    store_ACR_blast_region_in_syntenic_dic = {}
    count = 0
    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count == 1:
                final_line = eachline + '\t' + spe1_prefix + '_CelltypeCate' + '\t' +  'Syntenic_regionID'
                store_final_line_list.append(final_line)
            else:
                spe1ACR = col[1]
                store_ACR_blast_region_in_syntenic_dic[spe1ACR] = 1

                if spe1ACR in store_spe1_ACR_celltypecate_dic:
                    celltype_cate = store_spe1_ACR_celltypecate_dic[spe1ACR]
                else:
                    celltype_cate = 'Wrong'

                if spe1ACR in store_spe1_ACR_syntenic_region_dic:
                    syntenic_regionID_str = ','.join(store_spe1_ACR_syntenic_region_dic[spe1ACR])
                else:
                    syntenic_regionID_str = 'NotMoreThanTwoGenesInPair'

                final_line = eachline + '\t' + celltype_cate + '\t' + syntenic_regionID_str
                store_final_line_list.append(final_line)

    for eachACR in store_spe1_ACR_syntenic_region_dic:
        if eachACR not in store_ACR_blast_region_in_syntenic_dic:

            if eachACR in store_spe1_acr_spe1_ID_dic:
                ACRID = store_spe1_acr_spe1_ID_dic[eachACR]
            else:
                ACRID = 'none'

            if eachACR in store_spe1_ACR_celltypecate_dic:
                celltype_cate = store_spe1_ACR_celltypecate_dic[eachACR]
            else:
                celltype_cate = 'Wrong'

            if eachACR in store_target_ACR_cate_spe1_dic:
                InorNotH3K27m3_spe1 = store_target_ACR_cate_spe1_dic[eachACR]['InorNotH3K27m3']
                BroOrRes_spe1 = store_target_ACR_cate_spe1_dic[eachACR]['BroOrRes']
            else:
                InorNotH3K27m3_spe1 = 'none'
                BroOrRes_spe1 = 'none'

            syntenic_regionID_str = ','.join(store_spe1_ACR_syntenic_region_dic[eachACR])

            final_line = ACRID + '\t' + eachACR + '\t' + InorNotH3K27m3_spe1 + '\t' + BroOrRes_spe1 + '\t' + \
                         'none' + '\t' + 'none' + '\t' + 'none' + '\t' + 'none' + '\t' + celltype_cate + '\t' + syntenic_regionID_str
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




##########
##updating 102623 we will check if the gene pair or syntenic region having the same gene order and label to the final bed file
def subfunction_check_syntenic_region_gene_pair_direction (input_spe1_gff_fl,input_spe2_gff_fl,
                                                           input_spe1_syntenic_genes_all_os_ACRs_fl,input_spe2_syntenic_genes_all_os_ACRs_fl,
                                                           input_spe2_syntenic_small_region_fl,
                                                           input_summary_fl,opt_dir,input_spe2_acr_fl,
                                                           input_spe1_acr_celltype_decidedbyCover_fl,input_spe2_acr_celltype_decidebyCover_fl,
                                                           spe1_prefix,spe2_prefix):

    ##store the gene information
    store_gene_st_ed_dir_dic = {}
    with open (input_spe1_gff_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split()

                if col[2] == 'gene':

                    if spe1_prefix == 'rice':

                        if re.match('Name=(.+)',col[8].split(';')[1]):
                            mt = re.match('Name=(.+)',col[8].split(';')[1])
                            geneID = mt.group(1)
                            chrnm = col[0]
                            st = col[3]
                            ed = col[4]
                            dir = col[6]
                            store_gene_st_ed_dir_dic[geneID] = {'chrnm':chrnm,'st':st,'ed':ed,'dir':dir}

                    else:

                        if spe1_prefix == 'sorghum':
                            print(col[8].split(';')[0])
                            if re.match('ID=.+_pg\d+\..+',col[8].split(';')[0]):
                                mt = re.match('ID=(.+_pg\d+)\..+',col[8].split(';')[0])
                                geneID = mt.group(1)
                            else:

                                mt = re.match('ID=(.+)', col[8].split(';')[0])
                                geneID = mt.group(1)

                            chrnm = col[0]
                            st = col[3]
                            ed = col[4]
                            dir = col[6]
                            store_gene_st_ed_dir_dic[geneID] = {'chrnm': chrnm, 'st': st, 'ed': ed, 'dir': dir}

                        else:
                            if spe1_prefix == 'Uf':
                                if re.match('Name=(.+)',col[8].split(';')[1]):
                                    mt = re.match('Name=(.+)',col[8].split(';')[1])
                                    geneID = mt.group(1)
                                    chrnm = col[0]
                                    st = col[3]
                                    ed = col[4]
                                    dir = col[6]
                                    store_gene_st_ed_dir_dic[geneID] = {'chrnm': chrnm, 'st': st, 'ed': ed, 'dir': dir}

                            else:
                                if re.match('ID=(.+)',col[8].split(';')[0]):
                                    mt = re.match('ID=(.+)',col[8].split(';')[0])
                                    geneID = mt.group(1)
                                    chrnm = col[0]
                                    st = col[3]
                                    ed = col[4]
                                    dir = col[6]
                                    store_gene_st_ed_dir_dic[geneID] = {'chrnm':chrnm,'st':st,'ed':ed,'dir':dir}




    store_spe2_gene_st_ed_dir_dic = {}
    with open (input_spe2_gff_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split()

                if col[2] == 'gene':

                    if spe2_prefix == 'sorghum':

                        print(col[8].split(';')[0])
                        if re.match('ID=.+_pg\d+\..+',col[8].split(';')[0]):
                            mt = re.match('ID=(.+_pg\d+)\..+',col[8].split(';')[0])
                            geneID = mt.group(1)
                        else:

                            mt = re.match('ID=(.+)', col[8].split(';')[0])
                            geneID = mt.group(1)

                        chrnm = col[0]
                        st = col[3]
                        ed = col[4]
                        dir = col[6]
                        store_spe2_gene_st_ed_dir_dic[geneID] = {'chrnm': chrnm, 'st': st, 'ed': ed, 'dir': dir}

                    else:

                        if spe2_prefix == 'Uf':
                            if re.match('Name=(.+)',col[8].split(';')[1]):
                                mt = re.match('Name=(.+)',col[8].split(';')[1])
                                geneID = mt.group(1)
                                chrnm = col[0]
                                st = col[3]
                                ed = col[4]
                                dir = col[6]
                                store_spe2_gene_st_ed_dir_dic[geneID] = {'chrnm': chrnm, 'st': st, 'ed': ed, 'dir': dir}


                        else:

                            if spe2_prefix == 'rice':

                                if re.match('Name=(.+)', col[8].split(';')[1]):
                                    mt = re.match('Name=(.+)', col[8].split(';')[1])
                                    geneID = mt.group(1)
                                    chrnm = col[0]
                                    st = col[3]
                                    ed = col[4]
                                    dir = col[6]
                                    store_spe2_gene_st_ed_dir_dic[geneID] = {'chrnm': chrnm, 'st': st, 'ed': ed, 'dir': dir}

                            else:
                                if re.match('ID=(.+)',col[8].split(';')[0]):
                                    mt = re.match('ID=(.+)',col[8].split(';')[0])
                                    geneID = mt.group(1)
                                    chrnm = col[0]
                                    st = col[3]
                                    ed = col[4]
                                    dir = col[6]
                                    store_spe2_gene_st_ed_dir_dic[geneID] = {'chrnm':chrnm,'st':st,'ed':ed,'dir':dir}





    #print(store_gene_st_ed_dir_dic)

    store_syntenic_regionID_gene_dic = {}
    with open (input_spe1_syntenic_genes_all_os_ACRs_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            syntenic_region_ID = col[4]
            geneID = col[3]
            genedir = store_gene_st_ed_dir_dic[geneID]['dir']

            if syntenic_region_ID in store_syntenic_regionID_gene_dic:
                store_syntenic_regionID_gene_dic[syntenic_region_ID][geneID] = {'chr':col[0],'st':col[1],'ed':col[2],'dir':genedir}
            else:
                store_syntenic_regionID_gene_dic[syntenic_region_ID] = {}
                store_syntenic_regionID_gene_dic[syntenic_region_ID][geneID] = {'chr': col[0], 'st': col[1],
                                                                                'ed': col[2], 'dir': genedir}

    store_syntenic_regionID_feature_dic = {}
    for eachregion in store_syntenic_regionID_gene_dic:
        geneID_dic = store_syntenic_regionID_gene_dic[eachregion]

        store_dir_dic = {}
        for eachgeneID in geneID_dic:
            genedir = geneID_dic[eachgeneID]['dir']
            store_dir_dic[genedir] = 1

        if len(list(store_dir_dic.keys())) == 1:
            gene_pair_dir = 'OneDirection'
        else:
            gene_pair_dir = 'TwoDirection'

        store_syntenic_regionID_feature_dic[eachregion] = {'geneID_dic':geneID_dic,'genePair_dir':gene_pair_dir}

    store_spe2_syntenicID_blastsamllRegion_dic = {}
    with open (input_spe2_syntenic_small_region_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            small_region = col[0] + '_' + col[1] + '_' + col[2]
            syntenic_ID = col[3].split(';')[-1]
            if syntenic_ID in store_spe2_syntenicID_blastsamllRegion_dic:
                store_spe2_syntenicID_blastsamllRegion_dic[syntenic_ID][small_region] = 1
            else:
                store_spe2_syntenicID_blastsamllRegion_dic[syntenic_ID] = {}
                store_spe2_syntenicID_blastsamllRegion_dic[syntenic_ID][small_region] = 1


    store_final_line_list = []
    count = 0
    with open (input_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:

                if col[9] != 'NotMoreThanTwoGenesInPair':

                    syntenic_region_list = col[9].split(',')
                    genePair_dir_list = []
                    keep_syntenic_reigon_list = []
                    for eachsyntenic_region in syntenic_region_list:

                        genePair_dir = store_syntenic_regionID_feature_dic[eachsyntenic_region]['genePair_dir']
                        genePair_dir_list.append(genePair_dir)

                        ##updating 103123 we will ignore the one direction and focus on the both of them
                        #if genePair_dir == 'OneDirection':
                        keep_syntenic_reigon_list.append(eachsyntenic_region)


                    if keep_syntenic_reigon_list == []:
                        final_update_syntenic_region = 'none'
                    else:
                        final_update_syntenic_region = ','.join(keep_syntenic_reigon_list)

                    final_line = eachline + '\t' + ','.join(genePair_dir_list) + '\t' + final_update_syntenic_region
                    store_final_line_list.append(final_line)

                else:
                    store_final_line_list.append(eachline + '\t' + 'none' + '\t' + 'none')

            else:
                final_line = eachline + '\t' + 'riceGene_Direction' + '\t' + 'Syntenic_regionID_update'
                store_final_line_list.append(final_line)

    with open(opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    store_final_line_list = []
    count = 0
    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            spe2_small_region = col[4]

            if count == 1:
                store_final_line_list.append(eachline + '\t' + 'IfRightSyntenic')
            else:

                if col[11] != 'none':
                    syntenic_region_list = col[11].split(',')

                    store_syntenitc_smallspe2_region_dic = {}
                    for eachregion in syntenic_region_list:

                        ##Here we would see some syntenic regions do not show in the   store_spe2_syntenicID_blastsamllRegion_dic
                        if eachregion in store_spe2_syntenicID_blastsamllRegion_dic:

                            small_region_dic = store_spe2_syntenicID_blastsamllRegion_dic[eachregion]
                            for eachsmallregion in small_region_dic:
                                store_syntenitc_smallspe2_region_dic[eachsmallregion] = 1

                    ##do the filration
                    if spe2_small_region in store_syntenitc_smallspe2_region_dic:
                        final_line = eachline + '\t' + 'Yes'
                    else:
                        final_line = eachline + '\t' + 'No'

                else:
                    final_line = eachline + '\t' + 'none'

                store_final_line_list.append(final_line)


    with open(opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##updating 102623
    ##we will check the syntenic region of spe2
    store_syntenicregion_chr_genes_direction_dic = {}
    with open (input_spe2_syntenic_genes_all_os_ACRs_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            syntenic_ID = col[4]
            geneID = col[3]
            chr = col[0]
            gene_dir = store_spe2_gene_st_ed_dir_dic[geneID]['dir']

            if syntenic_ID in store_syntenicregion_chr_genes_direction_dic:
                if chr in store_syntenicregion_chr_genes_direction_dic[syntenic_ID]:
                    store_syntenicregion_chr_genes_direction_dic[syntenic_ID][chr][geneID] = {'dir':gene_dir}
                else:
                    store_syntenicregion_chr_genes_direction_dic[syntenic_ID][chr] = {}
                    store_syntenicregion_chr_genes_direction_dic[syntenic_ID][chr][geneID] = {'dir': gene_dir}
            else:
                store_syntenicregion_chr_genes_direction_dic[syntenic_ID] = {}
                store_syntenicregion_chr_genes_direction_dic[syntenic_ID][chr] = {}
                store_syntenicregion_chr_genes_direction_dic[syntenic_ID][chr][geneID] = {'dir': gene_dir}


    store_final_line_list = []
    count = 0
    with open(opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                spe2_small_region = col[4]
                #print(spe2_small_region)

                if spe2_small_region != 'none':
                    mt = re.match('(.+)_.+_.+',spe2_small_region)
                    spe2_chr = mt.group(1)

                    final_syntenic_str = col[11]
                    if final_syntenic_str == 'none':
                        final_line = eachline + '\t' + 'none'
                        store_final_line_list.append(final_line)
                    else:

                        store_final_cate_list = []
                        for eachsyntenic_region in final_syntenic_str.split(','):

                            final_dir_cate = 'Wrong'
                            if eachsyntenic_region in store_syntenicregion_chr_genes_direction_dic:


                                for eachchr in store_syntenicregion_chr_genes_direction_dic[eachsyntenic_region]:

                                    geneID_dir_dic = store_syntenicregion_chr_genes_direction_dic[eachsyntenic_region][eachchr]

                                    if eachchr == spe2_chr:
                                        ##we allow the target chrome has two genes
                                        if len(list(geneID_dir_dic.keys())) == 2:

                                            store_dir_dic = {}
                                            for eachgeneID in geneID_dir_dic:
                                                dir = geneID_dir_dic[eachgeneID]['dir']
                                                store_dir_dic[dir] = 1

                                            if len(list(store_dir_dic.keys())) == 1:
                                                final_dir_cate = 'OneDirection'
                                            else:
                                                final_dir_cate = 'TwoDirection'

                                        else:
                                            final_dir_cate = 'SinlgeGene'

                            else:
                                final_dir_cate = 'none'

                            store_final_cate_list.append(final_dir_cate)

                        store_final_cate_str = ','.join(store_final_cate_list)

                        final_line = eachline + '\t' + store_final_cate_str
                        store_final_line_list.append(final_line)

                else:
                    final_line = eachline + '\t' + 'none'
                    store_final_line_list.append(final_line)
            else:
                final_line = eachline + '\t' + spe2_prefix + 'Gene_Direction'
                store_final_line_list.append(final_line)


    with open(opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##upadting 110123
    ##add the cell type for the maize
    store_spe2_ACR_celltypecate_dic = {}
    with open (input_spe2_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split ()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            mt = re.match('(.+);(.+)',col[3])
            acrID = mt.group(1)
            celltypecate = mt.group(2)
            store_spe2_ACR_celltypecate_dic[acrnm] = celltypecate

    store_final_line_list = []
    count = 0
    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                spe2_acr = col[5]
                if spe2_acr in store_spe2_ACR_celltypecate_dic:
                    celltypecate = store_spe2_ACR_celltypecate_dic[spe2_acr]
                else:
                    celltypecate = 'none'

                final_line = eachline + '\t' + celltypecate
            else:
                final_line = eachline + '\t' + spe2_prefix + '_CelltypeCate'

            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##udpating 110223
    ##we will add the other two
    store_acr_spe1_celltypeCoverage_dic = {}
    with open (input_spe1_acr_celltype_decidedbyCover_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            celltype_str = col[4]
            celltype_cate = col[3]
            store_acr_spe1_celltypeCoverage_dic[acrnm] = {'acrnm':acrnm,'celltype_str':celltype_str,'celltype_cate':celltype_cate}

    store_acr_spe2_celltypeCoverage_dic = {}
    with open(input_spe2_acr_celltype_decidebyCover_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            celltype_str = col[4]
            celltype_cate = col[3]
            store_acr_spe2_celltypeCoverage_dic[acrnm] = {'acrnm': acrnm, 'celltype_str': celltype_str,
                                                     'celltype_cate': celltype_cate}

    store_final_line_list = []
    count = 0
    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                spe1_acr = col[1]
                spe2_acr = col[5]

                if spe1_acr in store_acr_spe1_celltypeCoverage_dic:
                    spe1_celltype_str = store_acr_spe1_celltypeCoverage_dic[spe1_acr]['celltype_str']
                    spe1_celltype_cate = store_acr_spe1_celltypeCoverage_dic[spe1_acr]['celltype_cate']
                else:
                    spe1_celltype_str = 'none'
                    spe1_celltype_cate = 'none'

                if spe2_acr in store_acr_spe2_celltypeCoverage_dic:
                    spe2_celltype_str = store_acr_spe2_celltypeCoverage_dic[spe2_acr]['celltype_str']
                    spe2_celltype_cate = store_acr_spe2_celltypeCoverage_dic[spe2_acr]['celltype_cate']
                else:
                    spe2_celltype_str = 'none'
                    spe2_celltype_cate = 'none'

                final_line = eachline + '\t' + spe1_celltype_str + '\t' + spe1_celltype_cate + '\t' + spe2_celltype_str + '\t' + spe2_celltype_cate
                store_final_line_list.append(final_line)

            else:
                first_line = eachline + '\t' + spe1_prefix + '_celltypeStr_CoverVer' + '\t' + spe1_prefix + '_celltypeCate_CoverVer' + '\t' + \
                             spe2_prefix + '_celltypeStr_CoverVer' + '\t' + spe2_prefix + '_celltypeCate_CoverVer'
                store_final_line_list.append(first_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


def subfunction_add_species_specific_one_under_non_syntenic_regions (ipt_final_summary_fl,spe1_H3K27me3_fl,spe1_acr_fl,
                                                                     spe1_prefix,spe2_prefix, opt_dir):

    store_all_acr_loc_dic = {}
    with open (spe1_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            store_all_acr_loc_dic[acrloc] = 1


    store_acr_loc_summary_dic = {}
    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                acr_loc = col[1]
                store_acr_loc_summary_dic[acr_loc] = 1

    store_acr_not_in_summary_dic = {}
    for eachacrloc in store_all_acr_loc_dic:
        if eachacrloc not in store_acr_loc_summary_dic:
            store_acr_not_in_summary_dic[eachacrloc] = 1

    ##for the H3K27me3 file
    store_target_ACR_cate_spe1_dic = {}
    count = 0
    with open(spe1_H3K27me3_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            broad_restrict_cate = col[3]
            final_cate = col[-1]
            count += 1
            if count != 1:
                store_target_ACR_cate_spe1_dic[acrnm] = {'BroOrRes': broad_restrict_cate,
                                                         'InorNotH3K27m3': final_cate}


    store_spe1_ACR_celltypecate_dic = {}
    with open(spe1_acr_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            mt = re.match('(.+);(.+)', col[3])
            acrID = mt.group(1)
            celltypecate = mt.group(2)
            store_spe1_ACR_celltypecate_dic[acrnm] = celltypecate


    store_final_line_list = []
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            store_final_line_list.append(eachline)

    for eachacr in store_acr_not_in_summary_dic:

        if eachacr in store_target_ACR_cate_spe1_dic:
            BroOrRes = store_target_ACR_cate_spe1_dic[eachacr]['BroOrRes']
            InorNotH3K27m3 = store_target_ACR_cate_spe1_dic[eachacr]['InorNotH3K27m3']
        else:
            BroOrRes = 'none'
            InorNotH3K27m3 = 'none'

        if eachacr in store_spe1_ACR_celltypecate_dic:
            celltypecate_spe1 = store_spe1_ACR_celltypecate_dic[eachacr]
        else:
            celltypecate_spe1 = 'none'

        none_list = []
        for i in range(10):
            none_list.append('none')

        final_line = 'none' + '\t' + eachacr + '\t' + InorNotH3K27m3 + '\t' + BroOrRes + '\t' + \
                     'none' + '\t' + 'none' + '\t' + 'none' + '\t' + 'none' + '\t' + \
                     celltypecate_spe1 + '\t' + '\t'.join(none_list)
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype_addSpeSpec.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



##updating 121123
def subfunction_check_not_shared_statistic_different_enrich_for_nonsyntenic_region (ipt_sharednotsharedCate_celltype_fl,opt_dir,spe1_prefix,spe2_prefix):

    ##we will check whether exact cell types were more enrich in the spec
    store_celltype_cate_acr_num_dic = {}
    with open (ipt_sharednotsharedCate_celltype_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acr = col[0]
            cate = col[2]
            celltype = col[1]

            if 'unknown_cells' not in celltype:

                if cate != 'NotBlast':

                    if celltype in store_celltype_cate_acr_num_dic:

                        if cate in store_celltype_cate_acr_num_dic[celltype]:
                            store_celltype_cate_acr_num_dic[celltype][cate][acr] = 1
                        else:
                            store_celltype_cate_acr_num_dic[celltype][cate] = {}
                            store_celltype_cate_acr_num_dic[celltype][cate][acr] = 1

                    else:
                        store_celltype_cate_acr_num_dic[celltype] = {}
                        store_celltype_cate_acr_num_dic[celltype][cate] = {}
                        store_celltype_cate_acr_num_dic[celltype][cate][acr] = 1

    ##for each cell type we conduct the enrichment for the target to anyother cell types
    store_final_line_list = []
    for eachcelltype in store_celltype_cate_acr_num_dic:

        ##for the target celltype enrichment
        ##we focuse on the not shared cate
        not_shared_Tcelltype_acr_num = len(list(store_celltype_cate_acr_num_dic[eachcelltype]['NotShared'].keys()))
        sharedA_Tcelltype_acr_num = len(list(store_celltype_cate_acr_num_dic[eachcelltype]['SharedAcc'].keys()))
        sharedI_Tcelltype_acr_num = len(list(store_celltype_cate_acr_num_dic[eachcelltype]['SharedInAcc'].keys()))

        Not_notshared_Tcelltype_acr_num = sharedA_Tcelltype_acr_num + sharedI_Tcelltype_acr_num

        for eachothercelltype in store_celltype_cate_acr_num_dic:



            not_shared_Ocelltype_acr_num = len(list(store_celltype_cate_acr_num_dic[eachothercelltype]['NotShared'].keys()))
            sharedA_Ocelltype_acr_num = len(list(store_celltype_cate_acr_num_dic[eachothercelltype]['SharedAcc'].keys()))
            sharedI_Ocelltype_acr_num = len(list(store_celltype_cate_acr_num_dic[eachothercelltype]['SharedInAcc'].keys()))

            Not_notshared_Ocelltype_acr_num = sharedA_Ocelltype_acr_num + sharedI_Ocelltype_acr_num

            ##do the fisher
            oddsratio, pvalue = stats.fisher_exact(
                [[not_shared_Tcelltype_acr_num, Not_notshared_Tcelltype_acr_num], [not_shared_Ocelltype_acr_num, Not_notshared_Ocelltype_acr_num]],
                alternative='greater')

            final_line = eachcelltype + '\t' + eachothercelltype + '\t' + str(pvalue) + '\t' + str(not_shared_Tcelltype_acr_num) + '\t' + \
                         str(Not_notshared_Tcelltype_acr_num) + '\t' + str(not_shared_Ocelltype_acr_num) + '\t' + str(Not_notshared_Ocelltype_acr_num)
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt2_enrichment_celltype_diff_cate_' + spe1_prefix + '_' + spe2_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


















