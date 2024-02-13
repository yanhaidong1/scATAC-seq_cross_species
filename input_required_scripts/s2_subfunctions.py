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


def subfunction_summarize_overview (ipt_final_summary_blast_fl,ipt_sep2_syntenic_region_fl,input_rice_phy_score_fl,opt_dir,
                                    ipt_spe1_ACR_gene_Cate_fl,ipt_spe2_ACR_gene_Cate_fl,ipt_spe1_ACR_celltype_fl,
                                    spe1_prefix,spe2_prefix,
                                    s2_s1_decide_celltype_cate,
                                    s2_s1_open_check_phy_score):


    store_syntenic_region_allRiceACR_dic = {}
    store_syntenic_region_RiceACR_MaizeACR_dic = {}
    store_allRiceACR_dic = {}
    store_allRiceACR_BlastToMaizeACR_dic = {}
    store_allRiceACR_BlastToMaizeRegionNotACR_dic = {}

    store_allMaizeACR_couldBeBlastedToRiceACR_dic = {}

    ##store CT acr and broad ACR
    store_allRiceACR_diffcate_dic = {}
    store_allRiceACR_BlastToMaizeACR_diffcate_dic = {}
    store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic = {}

    ##store maize acr
    store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic = {}

    ##store rice acr
    store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic = {}

    ##store rice and maize acr at same time
    ##udpating 111523
    store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic = {}


    ##store conserved broad and also the CT specific acr to check the score
    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic = {}

    ##store broad ACR to maize different cell type
    store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic = {}


    ##store_riceACR_ct_cate_dic
    ##updating 112523
    store_riceACR_CT_cate_dic = {}

    count = 0
    with open (ipt_final_summary_blast_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                riceACRID = col[0]
                riceACRloc = col[1]
                maizeACRloc = col[5]
                ricemaizeACRloc = riceACRloc + ',' + maizeACRloc

                acrline = '\t'.join(col[1].split('_'))
                store_riceACR_CT_cate_dic[acrline] = col[8]


                if s2_s1_decide_celltype_cate == 'CTVer':
                    rice_celltype = col[8]
                    maize_celltype = col[14]
                else:
                    if s2_s1_decide_celltype_cate == 'CoverVer':
                        rice_celltype = col[15]
                        maize_celltype = col[17]
                    else:
                        print('No rice and maize cell type identfied, please check the decide celltype cate is right or not')
                        break


                SyntenicRegionID_str = col[11]
                SyntenicRegionID_list = SyntenicRegionID_str.split(',')
                for eachSyntenicRegionID in SyntenicRegionID_list:
                    store_syntenic_region_allRiceACR_dic[eachSyntenicRegionID] = 1
                    store_allRiceACR_dic[riceACRloc] = 1

                    ##add cell type information
                    if s2_s1_decide_celltype_cate == 'CTVer':
                        if rice_celltype == 'broadly_accessible':
                            final_cate = 'broadACR'
                        else:
                            final_cate = 'restrictedACR'
                    else:

                        if s2_s1_decide_celltype_cate == 'CoverVer':
                            rice_celltype_cover_cate = col[16]
                            if rice_celltype_cover_cate == 'broad':
                                final_cate = 'broadACR'
                            else:
                                if rice_celltype_cover_cate == 'restricted':
                                    final_cate = 'restrictedACR'
                                else:
                                    final_cate = 'otherACR'

                        else:
                            print('please check the decide cell type cate is right or not')
                            break

                    if final_cate in store_allRiceACR_diffcate_dic:
                        store_allRiceACR_diffcate_dic[final_cate][riceACRloc] = 1
                    else:
                        store_allRiceACR_diffcate_dic[final_cate] = {}
                        store_allRiceACR_diffcate_dic[final_cate][riceACRloc] = 1


                    if riceACRID != 'none':

                        if maizeACRloc != 'none':

                            store_allRiceACR_BlastToMaizeACR_dic[riceACRloc] = 1
                            store_allMaizeACR_couldBeBlastedToRiceACR_dic[maizeACRloc] = 1

                            ##add cell type information
                            if s2_s1_decide_celltype_cate == 'CTVer':
                                if rice_celltype == 'broadly_accessible':
                                    final_rice_cate = 'broadACR'
                                else:
                                    final_rice_cate = 'restrictedACR'
                            else:
                                if s2_s1_decide_celltype_cate == 'CoverVer':
                                    rice_celltype_cover_cate = col[16]
                                    if rice_celltype_cover_cate == 'broad':
                                        final_rice_cate = 'broadACR'
                                    else:
                                        if rice_celltype_cover_cate == 'restricted':
                                            final_rice_cate = 'restrictedACR'
                                        else:
                                            final_rice_cate = 'otherACR'
                                else:
                                    print('please check the decide cell type cate is right or not')
                                    break


                            if final_rice_cate in store_allRiceACR_BlastToMaizeACR_diffcate_dic:
                                store_allRiceACR_BlastToMaizeACR_diffcate_dic[final_rice_cate][riceACRloc] = 1
                            else:
                                store_allRiceACR_BlastToMaizeACR_diffcate_dic[final_rice_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_diffcate_dic[final_rice_cate][riceACRloc] = 1

                            if s2_s1_decide_celltype_cate == 'CTVer':
                                if maize_celltype == 'broadly_accessible':
                                    final_maize_cate = 'broadACR'
                                else:
                                    final_maize_cate = 'restrictedACR'
                            else:
                                if s2_s1_decide_celltype_cate == 'CoverVer':
                                    maize_celltype_cover_cate = col[18]
                                    if maize_celltype_cover_cate == 'broad':
                                        final_maize_cate = 'broadACR'
                                    else:
                                        if maize_celltype_cover_cate == 'restricted':
                                            final_maize_cate = 'restrictedACR'
                                        else:
                                            final_maize_cate = 'otherACR'
                                else:
                                    print('please check the decide cell type cate is right or not')
                                    break


                            ##calculate the number of maize ACR in different cate
                            if final_rice_cate in store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic:
                                if final_maize_cate in store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate]:
                                    store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate][final_maize_cate][maizeACRloc] = 1
                                else:
                                    store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate][final_maize_cate] = {}
                                    store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate][final_maize_cate][maizeACRloc] = 1
                            else:
                                store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate][
                                    final_maize_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[final_rice_cate][
                                    final_maize_cate][maizeACRloc] = 1

                            ##calculate the number of rice ACR in different cate
                            if final_rice_cate in store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic:
                                if final_maize_cate in store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate]:
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate][final_maize_cate][riceACRloc] = 1
                                else:
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate][final_maize_cate] = {}
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate][final_maize_cate][riceACRloc] = 1
                            else:
                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate][final_maize_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[final_rice_cate][final_maize_cate][riceACRloc] = 1

                            ##updating 111523
                            ##store teh rice maize ACR in different cate
                            if final_rice_cate in store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic:
                                if final_maize_cate in store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[
                                    final_rice_cate]:
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate][
                                        final_maize_cate][ricemaizeACRloc] = 1
                                else:
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate][
                                        final_maize_cate] = {}
                                    store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate][
                                        final_maize_cate][ricemaizeACRloc] = 1
                            else:
                                store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate][final_maize_cate] = {}
                                store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[final_rice_cate][final_maize_cate][ricemaizeACRloc] = 1


                            ##updating 110223
                            ##we will check the ACR in different CT
                            if s2_s1_decide_celltype_cate == 'CTVer':
                                if rice_celltype != 'broadly_accessible' and maize_celltype != 'broadly_accessible' and \
                                        maize_celltype != 'procambial_meristem' and maize_celltype != 'procambium':

                                    rice_celltype_list = rice_celltype.split(',')
                                    maize_celltype_list = maize_celltype.split(',')

                                    for eachrice_celltype in rice_celltype_list:
                                        for eachmaize_celltype in maize_celltype_list:

                                            if eachrice_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic:
                                                if eachmaize_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachrice_celltype]:
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachrice_celltype][eachmaize_celltype][riceACRloc] = 1
                                                else:
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                        eachrice_celltype][eachmaize_celltype] = {}
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                        eachrice_celltype][eachmaize_celltype][riceACRloc] = 1

                                            else:
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachrice_celltype] = {}
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachrice_celltype][eachmaize_celltype] = {}
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachrice_celltype][eachmaize_celltype][riceACRloc] = 1

                                    ##updating 120923
                                    ##we will update the way to check the rice ACR loc corresponding to two
                                    #for eachrice_celltype in rice_celltype_list:

                                    #    if eachrice_celltype == 'companion_cell':
                                    #        rice_new_celltype = 'companion_cells_sieve_elements'
                                    #    else:
                                    #        rice_new_celltype = eachrice_celltype

                                    #    if rice_new_celltype in maize_celltype_list:
                                    #        rice_t_celltype = rice_new_celltype
                                    #    else:
                                    #        rice_t_celltype = rice_new_celltype + '__toOthers'

                                    #    if rice_t_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic:
                                    #        store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                    #            rice_t_celltype][riceACRloc] = 1
                                    #    else:
                                    #        store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                    #            rice_t_celltype] = {}
                                    #        store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                    #            rice_t_celltype][riceACRloc] = 1



                                ##updating 110223
                                if rice_celltype == 'broadly_accessible':

                                    maize_celltype_list = maize_celltype.split(',')

                                    for eachmaize_celltype in maize_celltype_list:

                                        if rice_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic:
                                            if eachmaize_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[rice_celltype]:
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                    rice_celltype][eachmaize_celltype][riceACRloc] = 1
                                            else:
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                    rice_celltype][eachmaize_celltype] = {}
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                    rice_celltype][eachmaize_celltype][riceACRloc] = 1
                                        else:
                                            store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                rice_celltype] = {}
                                            store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                rice_celltype][eachmaize_celltype] = {}
                                            store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[
                                                rice_celltype][eachmaize_celltype][riceACRloc] = 1


                            if s2_s1_decide_celltype_cate == 'CoverVer':

                                rice_celltype_cover_cate = col[16]
                                if rice_celltype_cover_cate == 'restricted':

                                    rice_celltype_list = rice_celltype.split(',')
                                    maize_celltype_list = maize_celltype.split(',')

                                    for eachrice_celltype in rice_celltype_list:
                                        for eachmaize_celltype in maize_celltype_list:

                                            if eachrice_celltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic:
                                                if eachmaize_celltype in \
                                                        store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                            eachrice_celltype]:
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                        eachrice_celltype][eachmaize_celltype][riceACRloc] = 1
                                                else:
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                        eachrice_celltype][eachmaize_celltype] = {}
                                                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                        eachrice_celltype][eachmaize_celltype][riceACRloc] = 1

                                            else:
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                    eachrice_celltype] = {}
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                    eachrice_celltype][eachmaize_celltype] = {}
                                                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[
                                                    eachrice_celltype][eachmaize_celltype][riceACRloc] = 1



                            ##for the syntenic check
                            if eachSyntenicRegionID in store_syntenic_region_RiceACR_MaizeACR_dic:

                                if riceACRID in store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID]:
                                    store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID][riceACRID][maizeACRloc] = 1
                                else:
                                    store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID][riceACRID] = {}
                                    store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID][riceACRID][maizeACRloc] = 1

                            else:
                                store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID] = {}
                                store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID][riceACRID] = {}
                                store_syntenic_region_RiceACR_MaizeACR_dic[eachSyntenicRegionID][riceACRID][maizeACRloc] = 1



                        else:
                            store_allRiceACR_BlastToMaizeRegionNotACR_dic[riceACRloc] = 1

                            ##add cell type information
                            if s2_s1_decide_celltype_cate == 'CTVer':
                                if rice_celltype == 'broadly_accessible':
                                    final_cate = 'broadACR'
                                else:
                                    final_cate = 'restrictedACR'
                            else:

                                if s2_s1_decide_celltype_cate == 'CoverVer':
                                    rice_celltype_cover_cate = col[16]
                                    if rice_celltype_cover_cate == 'broad':
                                        final_cate = 'broadACR'
                                    else:
                                        if rice_celltype_cover_cate == 'restricted':
                                            final_cate = 'restrictedACR'
                                        else:
                                            final_cate = 'otherACR'


                            if final_cate in store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic:
                                store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic[final_cate][riceACRloc] = 1
                            else:
                                store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic[final_cate] = {}
                                store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic[final_cate][riceACRloc] = 1


    store_total_syntenic_region_dic = {}
    with open (ipt_sep2_syntenic_region_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_total_syntenic_region_dic[col[3]] = 1

    syntenic_region_num = len(list(store_total_syntenic_region_dic.keys()))

    ##Here we will check if the syntenic region contains the acr
    syntenic_region_capturing_riceACR_dic = {}
    for eachregion in store_total_syntenic_region_dic:
        if eachregion in store_syntenic_region_allRiceACR_dic:
            syntenic_region_capturing_riceACR_dic[eachregion] = 1
    syntenic_region_real_num_capturing_riceACR = len(list(syntenic_region_capturing_riceACR_dic.keys()))

    ##check the diff between the real capture rice ACR
    for eachregion in store_syntenic_region_allRiceACR_dic:
        if eachregion not in syntenic_region_capturing_riceACR_dic:
            print(eachregion + ' has something wrong')

    ######
    ##opt1: check number of ACR and region for different case
    store_final_line_list = []
    syntenic_region_num_capturing_riceACR = len(list(store_syntenic_region_allRiceACR_dic.keys()))
    riceACR_num_within_syntenic_region = len(list(store_allRiceACR_dic.keys()))
    riceACR_num_within_syntenic_region_blastTomaizeACR = len(list(store_allRiceACR_BlastToMaizeACR_dic.keys()))
    riceACR_num_within_syntenic_region_blastTomaizeregionNotACR = len(list(store_allRiceACR_BlastToMaizeRegionNotACR_dic.keys()))

    final_line = 'syntenic_region_num' + '\t' + str(syntenic_region_num)
    store_final_line_list.append(final_line)
    final_line = 'syntenic_region_num_capturing_riceACR' + '\t' + str(syntenic_region_num_capturing_riceACR)
    store_final_line_list.append(final_line)
    final_line = 'syntenic_region_real_num_capturing_riceACR' + '\t' + str(syntenic_region_real_num_capturing_riceACR)
    store_final_line_list.append(final_line)
    final_line = 'riceACR_num_within_syntenic_region' + '\t' + str(riceACR_num_within_syntenic_region)
    store_final_line_list.append(final_line)
    final_line = 'riceACR_num_within_syntenic_region_blastTomaizeACR' + '\t' + str(riceACR_num_within_syntenic_region_blastTomaizeACR)
    store_final_line_list.append(final_line)
    final_line = 'riceACR_num_within_syntenic_region_blastTomaizeregionNotACR' + '\t' + str(riceACR_num_within_syntenic_region_blastTomaizeregionNotACR)
    store_final_line_list.append(final_line)

    ##we will store the cell type specific ones
    for eachcate in store_allRiceACR_diffcate_dic:
        final_line = eachcate + '_riceACR_num_within_syntenic_region' + '\t' + str(len(list(store_allRiceACR_diffcate_dic[eachcate].keys())))
        store_final_line_list.append(final_line)
    for eachcate in store_allRiceACR_BlastToMaizeACR_diffcate_dic:
        final_line = eachcate + '_riceACR_num_within_syntenic_region_blastTomaizeACR' + '\t' + str(len(list(store_allRiceACR_BlastToMaizeACR_diffcate_dic[eachcate].keys())))
        store_final_line_list.append(final_line)
    for eachcate in store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic:
        final_line = eachcate + '_riceACR_num_within_syntenic_region_blastTomaizeregionNotACR' + '\t' + str(len(list(store_allRiceACR_BlastToMaizeRegionNotACR_diffcate_dic[eachcate].keys())))
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_ACRnum.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##updating 111723 we will make a version to store all the ACRs
    store_final_acr_three_cate_line = []
    store_final_acr_three_cate_addCTcate_line = []
    for eachacr in store_allRiceACR_dic:
        acr_line = '\t'.join(eachacr.split('_'))
        final_line = acr_line + '\t' + 'All'
        store_final_acr_three_cate_line.append(final_line)

    for eachacr in store_allRiceACR_BlastToMaizeACR_dic:
        acr_line = '\t'.join(eachacr.split('_'))
        final_line =  acr_line + '\t' + 'SharedAcc'
        store_final_acr_three_cate_line.append(final_line)

    for eachacr in store_allRiceACR_BlastToMaizeRegionNotACR_dic:
        acr_line = '\t'.join(eachacr.split('_'))
        final_line = acr_line + '\t' + 'SharedInAcc'
        store_final_acr_three_cate_line.append(final_line)

    for eachacr in store_allRiceACR_dic:
        if eachacr not in store_allRiceACR_BlastToMaizeACR_dic:
            if eachacr not in store_allRiceACR_BlastToMaizeRegionNotACR_dic:
                acr_line = '\t'.join(eachacr.split('_'))
                final_line = acr_line + '\t' + 'NotShared'
                store_final_acr_three_cate_line.append(final_line)

    with open (opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs.txt','w+') as opt:
        for eachline in store_final_acr_three_cate_line:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs.txt > ' + \
          opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##store_riceACR_CT_cate_dic
    store_final_line_list = []
    with open (opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acr_line = col[0] + '\t' + col[1] + '\t' + col[2]
            CTcate = store_riceACR_CT_cate_dic[acr_line]
            final_line = eachline + '\t' + CTcate
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs_addCTcate.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs_addCTcate.txt > ' + \
          opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs_addCTcate_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##updating 111723
    ##we will generate an opt to store the cell type information
    ##this step will add the three category to the final acr file
    store_final_line_list = []
    with open (ipt_spe1_ACR_celltype_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]

            if acrnm in store_allRiceACR_dic:
                if acrnm not in store_allRiceACR_BlastToMaizeACR_dic:
                    if acrnm not in store_allRiceACR_BlastToMaizeRegionNotACR_dic:
                        final_cate = 'NotShared'
                    else:
                        final_cate = 'SharedInAcc'

                else:
                    final_cate = 'SharedAcc'

            else:
                final_cate = 'NotBlast'


            mt = re.match('.+;(.+)',col[3])
            celltype = mt.group(1)

            celltype_list = celltype.split(',')
            for eachcelltype in celltype_list:

                final_line = col[0] + '_' + col[1] + '_' + col[2] + '\t' + eachcelltype + '\t' + final_cate
                store_final_line_list.append(final_line)

    with open ( opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_celltypes.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



    ##updating 111523
    ##Here we will check if the broad syntenic broad were significantly enriched compared to the others
    ##use the fisher
    ##ACRtoACR_broad,  ACRtoACR_restricted
    ##ACRtoOthers_broad, ACRtoOthers_restricted

    ##oddsratio, pvalue = stats.fisher_exact([[8, 2], [1, 5]])

    ACRtoACR_broad = len(list(store_allRiceACR_BlastToMaizeACR_diffcate_dic['broadACR'].keys()))
    ACRtoACR_restricted = len(list(store_allRiceACR_BlastToMaizeACR_diffcate_dic['restrictedACR'].keys()))

    ACRtoOthers_broad = len(list(store_allRiceACR_diffcate_dic['broadACR'].keys())) - ACRtoACR_broad
    ACRtoOthers_restricted = len(list(store_allRiceACR_diffcate_dic['restrictedACR'].keys())) - ACRtoACR_restricted

    oddsratio, pvalue = stats.fisher_exact([[ACRtoACR_broad, ACRtoACR_restricted], [ACRtoOthers_broad, ACRtoOthers_restricted]],alternative='greater')

    with open (opt_dir + '/opt1_syntenic_region_' + spe1_prefix + '_ACRnum_check_ACRtoACR_enrich.txt','w+') as opt:
        final_line = 'pval' + '\t' + 'ACRtoACR_broad' + '\t' + 'ACRtoACR_restricted' + '\t' + 'ACRtoOthers_broad' + '\t' + 'ACRtoOthers_restricted'
        opt.write(final_line + '\n')
        final_line = str(pvalue) + '\t' + str(ACRtoACR_broad) + '\t' + str(ACRtoACR_restricted) + '\t' + str(ACRtoOthers_broad) + '\t' + str(ACRtoOthers_restricted)
        opt.write(final_line + '\n')








    ######
    ##opt2: check cate of ACR in maize to the cate of ACRs in rice
    ##updating 110123
    ##check how many restricted ACRs in rice would be the broad or restricted in maize
    ##updating 111523 we will add the category of close to genes
    store_spe1_acr_gene_cate_dic = {}
    with open (ipt_spe1_ACR_gene_Cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0]
            cate = col[1]
            store_spe1_acr_gene_cate_dic[acrloc] = cate

    store_spe2_acr_gene_cate_dic = {}
    with open (ipt_spe2_ACR_gene_Cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0]
            cate = col[1]
            store_spe2_acr_gene_cate_dic[acrloc] = cate


    store_final_line_list = []
    store_final_line_acr_ver_list = []
    first_line = spe1_prefix + 'CTcate' + '\t' + spe2_prefix + 'CTcate' + '\t' + spe2_prefix + 'ACRnum'
    store_final_line_list.append(first_line)
    for eachriceCTcate in store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic:
        for eachmaizeCTcate in store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[eachriceCTcate]:
            maize_acr_dic = store_allRiceACR_BlastToMaizeACR_TargetmaizeACR_diffcate_dic[eachriceCTcate][eachmaizeCTcate]
            final_line = eachriceCTcate + '\t' + eachmaizeCTcate + '\t' + str(len(list(maize_acr_dic.keys())))
            store_final_line_list.append(final_line)
            for eachmaizeacr in maize_acr_dic:
                final_line = '\t'.join(eachmaizeacr.split('_')) + '\t' + eachmaizeCTcate + '\t' + eachmaizeacr
                store_final_line_acr_ver_list.append(final_line)

    with open (opt_dir + '/opt2_syntenic_region_' + spe2_prefix + '_in_' + spe1_prefix + '_ACRnum_celltype.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')
    with open (opt_dir + '/opt2_syntenic_region_' + spe2_prefix + '_in_' + spe1_prefix + '_ACRName_celltype.txt','w+') as opt:
        for eachline in store_final_line_acr_ver_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt2_syntenic_region_' + spe2_prefix + '_in_' + spe1_prefix + '_ACRName_celltype.txt > ' + \
          opt_dir + '/opt2_syntenic_region_' + spe2_prefix + '_in_' + spe1_prefix + '_ACRName_celltype_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_final_line_list = []
    store_final_line_acr_ver_list = []
    first_line = spe1_prefix + 'CTcate' + '\t' + spe2_prefix + 'CTcate' + '\t' + spe1_prefix + 'ACRnum'
    store_final_line_list.append(first_line)
    for eachriceCTcate in store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic:
        for eachmaizeCTcate in store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[eachriceCTcate]:
            rice_acr_dic = store_allRiceACR_BlastToMaizeACR_TargetriceACR_diffcate_dic[eachriceCTcate][eachmaizeCTcate]
            final_line = eachriceCTcate + '\t' + eachmaizeCTcate + '\t' + str(len(list(rice_acr_dic.keys())))
            store_final_line_list.append(final_line)
            for eachriceacr in rice_acr_dic:
                if eachriceacr in store_spe1_acr_gene_cate_dic:
                    riceacr_gene_cate = store_spe1_acr_gene_cate_dic[eachriceacr]
                else:
                    riceacr_gene_cate = 'none'
                final_line = '\t'.join(eachriceacr.split('_')) + '\t' + eachriceCTcate + '\t' + eachmaizeCTcate + '\t' + riceacr_gene_cate
                store_final_line_acr_ver_list.append(final_line)

    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRnum_CelltypeCate.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')
    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate.txt','w+') as opt:
        for eachline in store_final_line_acr_ver_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate.txt > ' + \
          opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)


    ##udpating 111523
    ##opt2: check number and name of ACR pairs (for both Os and Zm)
    store_final_line_acrnum_list = []
    first_line = spe1_prefix + 'CTcate' + '\t' + spe2_prefix + 'CTcate' + '\t' + spe1_prefix + 'ACRnum' + '\t' + spe2_prefix + 'ACRnum'
    store_final_line_acrnum_list.append(first_line)

    store_final_line_acrloc_ver_list = []
    first_line = spe1_prefix + 'ACRchr' + '\t' + spe1_prefix + 'ACRst' + '\t' + spe1_prefix + 'ACRed' + '\t' + \
                 spe1_prefix + 'CTcate' + '\t' + spe1_prefix + 'GeneCate' + '\t' + \
                 spe2_prefix + 'ACRchr' + '\t' + spe2_prefix + 'ACRst' + '\t' + spe2_prefix + 'ACRed' + '\t' + \
                 spe2_prefix + 'CTcate' + '\t' + spe2_prefix + 'GeneCate'
    store_final_line_acrloc_ver_list.append(first_line)
    for eachriceCTcate in store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic:
        for eachmaizeCTcate in store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[eachriceCTcate]:

            ##store the number
            rice_maize_acr_dic = store_allRiceACR_BlastToMaizeACR_TargetriceACRmaizeACR_diffcate_dic[eachriceCTcate][eachmaizeCTcate]
            final_line = eachriceCTcate + '\t' + eachmaizeCTcate + '\t' + str(len(list(rice_maize_acr_dic.keys())))
            store_final_line_acrnum_list.append(final_line)

            ##store the loc
            for eachricemaizeacr in rice_maize_acr_dic:

                rice_acr = eachricemaizeacr.split(',')[0]
                maize_acr = eachricemaizeacr.split(',')[1]

                if rice_acr in store_spe1_acr_gene_cate_dic:
                    riceacr_gene_cate = store_spe1_acr_gene_cate_dic[rice_acr]
                else:
                    riceacr_gene_cate = 'none'

                if maize_acr in store_spe2_acr_gene_cate_dic:
                    maizeacr_gene_cate = store_spe2_acr_gene_cate_dic[maize_acr]
                else:
                    maizeacr_gene_cate = 'none'

                final_line = '\t'.join(rice_acr.split('_')) + '\t' + eachriceCTcate + '\t' + riceacr_gene_cate + '\t' + \
                             '\t'.join(maize_acr.split('_')) + '\t' + eachmaizeCTcate + '\t' + maizeacr_gene_cate
                store_final_line_acrloc_ver_list.append(final_line)

    with open(opt_dir + '/opt2_forBothACR_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRnum_CelltypeCate.txt','w+') as opt:
        for eachline in store_final_line_acrnum_list:
            opt.write(eachline + '\n')
    with open(opt_dir + '/opt2_forBothACR_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate.txt','w+') as opt:
        for eachline in store_final_line_acrloc_ver_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/opt2_forBothACR_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate.txt > ' + \
          opt_dir + '/opt2_forBothACR_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)



    ##updating 120423
    ##opt2: check number of shared-InAcc ACRs in different proximity to the nearby genes
    store_rice_CTtype_sharedIn_acrloc_dic = {}
    store_rice_CTtype_notshared_acrloc_dic = {}
    count = 0
    with open(ipt_final_summary_blast_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                riceACRID = col[0]
                riceACRloc = col[1]
                maizeACRloc = col[5]
                rice_celltype = col[8]

                if rice_celltype == 'broadly_accessible':
                    final_cate = 'Broad'
                else:
                    final_cate = 'CT'

                if riceACRID != 'none':
                    if maizeACRloc == 'none':

                        if final_cate in store_rice_CTtype_sharedIn_acrloc_dic:
                            store_rice_CTtype_sharedIn_acrloc_dic[final_cate][riceACRloc] = 1
                        else:
                            store_rice_CTtype_sharedIn_acrloc_dic[final_cate] = {}
                            store_rice_CTtype_sharedIn_acrloc_dic[final_cate][riceACRloc] = 1

                else:

                    if final_cate in store_rice_CTtype_notshared_acrloc_dic:
                        store_rice_CTtype_notshared_acrloc_dic[final_cate][riceACRloc] = 1
                    else:
                        store_rice_CTtype_notshared_acrloc_dic[final_cate] = {}
                        store_rice_CTtype_notshared_acrloc_dic[final_cate][riceACRloc] = 1



    store_final_line_list = []
    for eachfinalcate in store_rice_CTtype_sharedIn_acrloc_dic:
        riceacrloc_dic = store_rice_CTtype_sharedIn_acrloc_dic[eachfinalcate]

        for eachacrloc in riceacrloc_dic:

            rice_acr_proximity = store_spe1_acr_gene_cate_dic[eachacrloc]

            final_line = eachfinalcate + '\t' + eachacrloc + '\t' + rice_acr_proximity
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt2_for_ACR_sharedInAcc_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    store_final_line_list = []
    for eachfinalcate in store_rice_CTtype_notshared_acrloc_dic:
        riceacrloc_dic = store_rice_CTtype_notshared_acrloc_dic[eachfinalcate]

        for eachacrloc in riceacrloc_dic:

            rice_acr_proximity = store_spe1_acr_gene_cate_dic[eachacrloc]

            final_line = eachfinalcate + '\t' + eachacrloc + '\t' + rice_acr_proximity
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt2_for_ACR_notshared_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




    ##opt2: check cell type of ACRs in both rice and maize
    store_final_line_allacr_list = []
    store_final_line_acrcount_list = []
    for eachricecelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic:
        for eachmaizecelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachricecelltype]:
            riceacr_dic =  store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachricecelltype][eachmaizecelltype]
            riceacr_count = len(list(riceacr_dic.keys()))
            final_line = eachricecelltype + '\t' + eachmaizecelltype + '\t' + str(riceacr_count)
            store_final_line_acrcount_list.append(final_line)

            for eachacr in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachricecelltype][eachmaizecelltype]:
                final_line = '\t'.join(eachacr.split('_')) + '\t' + eachricecelltype + '\t' + eachmaizecelltype
                store_final_line_allacr_list.append(final_line)

    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_EachCelltype.txt','w+') as opt:
        for eachline in store_final_line_allacr_list:
            opt.write(eachline + '\n')

    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRnum_EachCelltype.txt','w+') as opt:
        for eachline in store_final_line_acrcount_list:
            opt.write(eachline + '\n')


    ##udpating 121023
    ##we will check teh name file and recaculate the name of ACR for the corresponding
    ##the idea is to check the same ACR with different spe pairs and decide which pair is not right
    store_acrloc_celltype_dic = {}
    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_EachCelltype.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            spe1_celltype = col[3]
            spe2_celltype = col[4]
            spe_celltype = spe1_celltype + '__' + spe2_celltype

            if acrloc in store_acrloc_celltype_dic:
                store_acrloc_celltype_dic[acrloc][spe_celltype] = 1
            else:
                store_acrloc_celltype_dic[acrloc] = {}
                store_acrloc_celltype_dic[acrloc][spe_celltype] = 1

    store_final_line_list = []
    for eachacrloc in store_acrloc_celltype_dic:

        celltype_pair_dic = store_acrloc_celltype_dic[eachacrloc]

        store_all_spe1_celltype_dic = {}
        store_all_spe2_celltype_dic = {}
        store_all_spe1_celltype_spe2_celltype_dic = {}
        for eachcelltypepair in celltype_pair_dic:
            mt = re.match('(.+)__(.+)',eachcelltypepair)
            spe1_celltype = mt.group(1)
            store_all_spe1_celltype_dic[spe1_celltype] = 1
            spe2_celltype = mt.group(2)
            store_all_spe2_celltype_dic[spe2_celltype] = 1

            if spe1_celltype in store_all_spe1_celltype_spe2_celltype_dic:
                store_all_spe1_celltype_spe2_celltype_dic[spe1_celltype][spe2_celltype] = 1
            else:
                store_all_spe1_celltype_spe2_celltype_dic[spe1_celltype] = {}
                store_all_spe1_celltype_spe2_celltype_dic[spe1_celltype][spe2_celltype] = 1


        for eachspe1_celltype in store_all_spe1_celltype_dic:

            if eachspe1_celltype in store_all_spe2_celltype_dic:
                check_spe1_celltype_have_correspd = 'yes'
            else:
                check_spe1_celltype_have_correspd = 'no'

            if check_spe1_celltype_have_correspd == 'yes':
                new_pair = eachspe1_celltype + '\t' + eachspe1_celltype
                final_line = '\t'.join(eachacrloc.split('_')) + '\t' + new_pair
                store_final_line_list.append(final_line)

            else:

                spe2_celltype_dic = store_all_spe1_celltype_spe2_celltype_dic[eachspe1_celltype]
                for eachspe2_celltype in spe2_celltype_dic:
                    new_pair = eachspe1_celltype + '\t' + eachspe2_celltype
                    final_line = '\t'.join(eachacrloc.split('_')) + '\t' + new_pair
                    store_final_line_list.append(final_line)


    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_EachCelltype_adjustPair.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##now calculate the new pair
    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic = {}
    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_EachCelltype_adjustPair.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            spe1_ct = col[3]
            spe2_ct = col[4]

            if spe1_ct in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic:
                if spe2_ct in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[spe1_ct]:
                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[spe1_ct][spe2_ct][acrloc] = 1
                else:
                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[spe1_ct][spe2_ct] = {}
                    store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[spe1_ct][spe2_ct][acrloc] = 1
            else:
                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[spe1_ct] = {}
                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[spe1_ct][spe2_ct] = {}
                store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[spe1_ct][spe2_ct][acrloc] = 1

    store_final_line_acrcount_list = []
    for eachricecelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic:
        for eachmaizecelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachricecelltype]:
            riceacr_dic = store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachricecelltype][eachmaizecelltype]
            riceacr_count = len(list(riceacr_dic.keys()))
            final_line = eachricecelltype + '\t' + eachmaizecelltype + '\t' + str(riceacr_count)
            store_final_line_acrcount_list.append(final_line)

    with open(opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRnum_EachCelltype_adjustPair.txt',
              'w+') as opt:
        for eachline in store_final_line_acrcount_list:
            opt.write(eachline + '\n')






    ##updating 120923
    ##does not work for the others
    ##check the cell type of ACRs in both rice and maize using the new way
    #store_final_line_allacr_list = []
    #store_final_line_acrcount_list = []

    #for eachriccelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic:

    #    riceacr_dic = store_allRiceACR_BlastToMaizeACR_TargetriceACR_CTcate_dic[eachriccelltype]
    #    riceacr_count = len(list(riceacr_dic.keys()))
    #    final_line = eachriccelltype + '\t' + str(riceacr_count)
    #    store_final_line_acrcount_list.append(final_line)

    #    for eachacr in riceacr_dic:
    #        final_line = '\t'.join(eachacr.split('_')) + '\t' + eachriccelltype
    #        store_final_line_allacr_list.append(final_line)

    #with open(opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_EachCelltype.txt',
    #          'w+') as opt:
    #    for eachline in store_final_line_allacr_list:
    #        opt.write(eachline + '\n')

    #with open(opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRnum_EachCelltype.txt',
    #          'w+') as opt:
    #    for eachline in store_final_line_acrcount_list:
    #        opt.write(eachline + '\n')







    ##updating 110223
    ##opt2: check the broad ACRs in both rice and maize
    store_final_line_allacr_list = []
    store_final_line_acrcount_list = []
    for eachricecelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic:
        for eachmaizecelltype in store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[eachricecelltype]:
            riceacr_dic = store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[eachricecelltype][eachmaizecelltype]
            riceacr_count = len(list(riceacr_dic.keys()))
            final_line = eachricecelltype + '\t' + eachmaizecelltype + '\t' + str(riceacr_count)
            store_final_line_acrcount_list.append(final_line)

            for eachacr in store_allRiceACR_BlastToMaizeACR_TargetriceACR_BroadToOthers_dic[eachricecelltype][eachmaizecelltype]:
                final_line = '\t'.join(eachacr.split('_')) + '\t' + eachricecelltype + '\t' + eachmaizecelltype
                store_final_line_allacr_list.append(final_line)

    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_BroadACRtoEachCelltype.txt','w+') as opt:
        for eachline in store_final_line_allacr_list:
            opt.write(eachline + '\n')

    with open (opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRnum_BroadACRtoEachCelltype.txt','w+') as opt:
        for eachline in store_final_line_acrcount_list:
            opt.write(eachline + '\n')


    #####################################################
    ##check the conesrvation score for these cate of ACRs
    if s2_s1_open_check_phy_score == 'yes':

        store_physcore_dir = opt_dir + '/store_physcore_dir'
        if not os.path.exists(store_physcore_dir):
            os.makedirs(store_physcore_dir)

        ipt_spe1_acr_cate_fl = opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_CelltypeCate_sorted.txt'

        cmd = 'bedtools intersect -wa -wb -a ' + ipt_spe1_acr_cate_fl + ' -b ' + input_rice_phy_score_fl + ' > ' + \
              store_physcore_dir + '/opt_intersect_' + spe1_prefix + '_acr_physcore_celltypecate.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ipt_spe1_acr_celltype_fl = opt_dir + '/opt2_syntenic_region_' + spe1_prefix + '_in_' + spe2_prefix + '_ACRName_EachCelltype.txt'

        cmd = 'bedtools intersect -wa -wb -a ' + ipt_spe1_acr_celltype_fl + ' -b ' + input_rice_phy_score_fl + ' > ' + \
              store_physcore_dir + '/opt_intersect_' + spe1_prefix + '_acr_physcore_percelltype.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##updating 110823
        ##compare the phylo score between CT and Broad peaks for alls
        store_rice_ct_broad_line_dic = {}
        count = 0
        with open(ipt_final_summary_blast_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                count += 1
                if count != 1:
                    riceACRID = col[0]
                    riceACRloc = col[1]
                    maizeACRloc = col[5]
                    ricecelltypecate = col[8]
                    maizecelltypecate = col[14]

                    if ricecelltypecate == 'broadly_accessible':
                        final_line = '\t'.join(riceACRloc.split('_')) + '\t' + ricecelltypecate
                        store_rice_ct_broad_line_dic[final_line] = 1
                    else:
                        final_line = '\t'.join(riceACRloc.split('_')) + '\t' + 'celltypeSpecific'
                        store_rice_ct_broad_line_dic[final_line] = 1

        with open (store_physcore_dir + '/temp_rice_broad_CT_acr.txt','w+') as opt:
            for eachline in store_rice_ct_broad_line_dic:
                opt.write(eachline + '\n')

        cmd = 'sort -k1,1V -k2,2n ' + store_physcore_dir + '/temp_rice_broad_CT_acr.txt > ' + \
              store_physcore_dir + '/temp_rice_broad_CT_acr_sorted.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

        cmd = 'bedtools intersect -wa -wb -a ' + store_physcore_dir + '/temp_rice_broad_CT_acr_sorted.txt' + ' -b ' + input_rice_phy_score_fl + ' > ' + \
              store_physcore_dir + '/opt_intersect_' + spe1_prefix + '_acr_physcore_broadCTACR.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)




    ######
    ##opt3: check one ACR in rice may correspond to how many ACRs in maize
    store_allRiceACR_blastToMaizeACR_num_dic = {}
    store_allRiceACR_blastToMaizeACR_line_list = []
    for eachsyntenic_region in store_syntenic_region_RiceACR_MaizeACR_dic:
        riceACRID_maizeACRloc_dic = store_syntenic_region_RiceACR_MaizeACR_dic[eachsyntenic_region]

        for eachriceACRID in riceACRID_maizeACRloc_dic:
            maizeACRloc_dic = riceACRID_maizeACRloc_dic[eachriceACRID]
            maizeACRnum = str(len(list(maizeACRloc_dic.keys())))
            maizeACRnum_name = 'BlastTo_' + maizeACRnum + '_' + spe2_prefix + 'ACRs'

            if maizeACRnum_name in store_allRiceACR_blastToMaizeACR_num_dic:
                store_allRiceACR_blastToMaizeACR_num_dic[maizeACRnum_name] += 1
            else:
                store_allRiceACR_blastToMaizeACR_num_dic[maizeACRnum_name] = 1

            maizeACRloc_str = ','.join(list(maizeACRloc_dic.keys()))

            num_maize_acr = len(list(maizeACRloc_dic.keys()))

            ##check if it is from the same or different col
            store_chrnm_dic = {}
            for eachmaizeacr in maizeACRloc_dic:
                mt = re.match('(.+)_.+_.+',eachmaizeacr)
                chrnm = mt.group(1)
                store_chrnm_dic[chrnm] = 1

            final_line = eachsyntenic_region + '\t' + eachriceACRID  + '\t' + str(num_maize_acr) + '\t' + str(len(store_chrnm_dic.keys())) + '\t' + maizeACRloc_str
            store_allRiceACR_blastToMaizeACR_line_list.append(final_line)


    with open (opt_dir + '/opt3_num_' + spe1_prefix + '_blastTo_diffnumOf_' + spe2_prefix + '.txt','w+') as opt:
        for eachcate in store_allRiceACR_blastToMaizeACR_num_dic:
            num_ACR = store_allRiceACR_blastToMaizeACR_num_dic[eachcate]
            final_line = eachcate + '\t' + str(num_ACR)
            opt.write(final_line + '\n')

    with open (opt_dir + '/opt3_' + spe1_prefix + '_blastTo_' + spe2_prefix + '.txt','w+') as opt:
        for eachline in store_allRiceACR_blastToMaizeACR_line_list:
            opt.write(eachline + '\n')


    ##updating 110323
    ##Here we will generate different files plotting to the browser
    store_regions_ACR_to_upload_To_browser_dir = opt_dir + '/store_regions_ACR_to_upload_To_browser_dir'
    if not os.path.exists(store_regions_ACR_to_upload_To_browser_dir):
        os.makedirs(store_regions_ACR_to_upload_To_browser_dir)

    #riceACR_num_within_syntenic_region = len(list(store_allRiceACR_dic.keys()))
    #riceACR_num_within_syntenic_region_blastTomaizeACR = len(list(store_allRiceACR_BlastToMaizeACR_dic.keys()))
    #riceACR_num_within_syntenic_region_blastTomaizeregionNotACR = len(list(store_allRiceACR_BlastToMaizeRegionNotACR_dic.keys()))

    ##for all the rice ACRs in the syntenic regions
    with open (store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion.txt','w+') as opt:
        for eachline in store_allRiceACR_dic:
            acr_line = '\t'.join(eachline.split('_'))
            opt.write(acr_line + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion.txt > ' + \
          store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##for ACR in syntenic region blastToMaize ACR
    with open (store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion_BlastToMaizeACR.txt','w+') as opt:
        for eachline in store_allRiceACR_BlastToMaizeACR_dic:
            acr_line = '\t'.join(eachline.split('_'))
            opt.write(acr_line + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion_BlastToMaizeACR.txt' + ' > ' + \
          store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe1_prefix + '_ACRs_in_' + spe1_prefix + '_SyntenicRegion_BlastToMaizeACR_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##we will check the spe2 ACR that could be blasted to rice ACR
    with open (store_regions_ACR_to_upload_To_browser_dir + '/opt_'+ spe2_prefix + '_ACRs_couldBeBlastedTo_' + spe1_prefix + '_ACRs.txt','w+') as opt:
        for eachline in store_allMaizeACR_couldBeBlastedToRiceACR_dic:
            acr_line = '\t'.join(eachline.split('_'))
            opt.write(acr_line + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + store_regions_ACR_to_upload_To_browser_dir + '/opt_'+ spe2_prefix + '_ACRs_couldBeBlastedTo_' + spe1_prefix + '_ACRs.txt > ' + \
          store_regions_ACR_to_upload_To_browser_dir + '/opt_' + spe2_prefix + '_ACRs_couldBeBlastedTo_' + spe1_prefix + '_ACRs_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)


##updating 121123
def subfunction_check_not_shared_statistic_different_enrich (ipt_sharednotsharedCate_celltype_fl,opt_dir,spe1_prefix,spe2_prefix):

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




##updating 111723
def subfunction_check_acc_in_rice_atlas_peaks (ipt_rice_atlas_peak_coverage_fl,ipt_total_rice_atlas_acr_fl,ipt_rice_leaf_syntenic_fl,ipt_rice_leaf_syntenic_addCTcate_fl,
                                               opt_dir,
                                               spe1_prefix,spe2_prefix,
                                               s2_s2_organct_coverage_cutoff):


    ##s2_s2_organct_coverage_cutoff would be 2.5 for the rice atlas

    store_all_rice_atlas_peak_dic = {}
    with open (ipt_total_rice_atlas_acr_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrnm = col[0] + '_' + col[1] + '_' + col[2]
            store_all_rice_atlas_peak_dic[acrnm] = 1


    ##we will check the H3K27m3 for each of organ one by one
    store_peak_organct_val_dic = {}
    store_total_organct_dic = {}
    with open(ipt_rice_atlas_peak_coverage_fl, 'r') as ipt:
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

    with open(opt_dir + '/temp_summary_DA_peak_string_filterCoverage.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##check the overlapping with the single rice called peaks
    ##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_cross_species_091823/output_dir_103123/step01_species_compare_add_cate_dir/store_rice_to_maize_dir
    ##Here we will use the three categories files
    cmd = 'bedtools intersect -wa -wb -a ' + ipt_rice_leaf_syntenic_fl + ' -b '  + ipt_total_rice_atlas_acr_fl + ' > ' + \
          opt_dir + '/temp_intersect_' + spe1_prefix + '_leafACR_to_atlas.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_final_line_list = []
    store_final_line_list_add_all_celltype_str = []
    with open (opt_dir + '/temp_intersect_' + spe1_prefix + '_leafACR_to_atlas.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            atlas_acr = col[4] + '_' + col[5] + '_' + col[6]

            celltypenum = store_eachpeak_celltype_str_num_dic[atlas_acr]['celltypenum']
            celltypestr = store_eachpeak_celltype_str_num_dic[atlas_acr]['celltypestr']

            final_line = eachline + '\t' + celltypenum
            store_final_line_list.append(final_line)

            final_line = eachline + '\t' + celltypenum + '\t' + celltypestr
            store_final_line_list_add_all_celltype_str.append(final_line)

    with open (opt_dir + '/opt_final_syntenic_ACR_' + spe1_prefix + '_' + spe2_prefix + '_addAtlasACRcelltype.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    with open (opt_dir + '/opt_final_syntenic_ACR_' + spe1_prefix + '_' + spe2_prefix + '_addAtlasACRcelltype_celltypestr.txt','w+') as opt:
        for eachline in store_final_line_list_add_all_celltype_str:
            opt.write(eachline + '\n')



    cmd = 'bedtools intersect -wa -wb -a ' + ipt_rice_leaf_syntenic_addCTcate_fl + ' -b '  + ipt_total_rice_atlas_acr_fl + ' > ' + \
          opt_dir + '/temp_intersect_' + spe1_prefix + '_leafACR_to_atlas_addCTcate.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_final_line_list = []
    store_final_line_addstr_list = []
    with open (opt_dir + '/temp_intersect_' + spe1_prefix + '_leafACR_to_atlas_addCTcate.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            atlas_acr = col[5] + '_' + col[6] + '_' + col[7]

            celltypenum = store_eachpeak_celltype_str_num_dic[atlas_acr]['celltypenum']

            final_line = eachline + '\t' + celltypenum
            store_final_line_list.append(final_line)

            celltypestr = store_eachpeak_celltype_str_num_dic[atlas_acr]['celltypestr']
            final_line = eachline + '\t' + celltypenum + '\t' + celltypestr
            store_final_line_addstr_list.append(final_line)

    with open (opt_dir + '/opt_final_syntenic_ACR_' + spe1_prefix + '_' + spe2_prefix + '_addAtlasACRcelltype_addCTcate.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    with open (opt_dir + '/opt_final_syntenic_ACR_' + spe1_prefix + '_' + spe2_prefix + '_addAtlasACRcelltype_CTstr_addCTcate.txt','w+') as opt:
        for eachline in store_final_line_addstr_list:
            opt.write(eachline + '\n')







##updating 120623
def subfunction_check_nearby_genes (ipt_final_summary_fl,ipt_spe1_gff_fl,ipt_spe1_acr_sumit_gene_closest_fl,ipt_spe1_gene_annot_fl,opt_dir,spe1_prefix,spe2_prefix):


    store_geneloc_geneID_dic = {}
    with open(ipt_spe1_gff_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split()
                if col[2] == 'gene':

                    chrnm = col[0]
                    st = col[3]
                    ed = col[4]

                    geneloc = chrnm + '_' + st + '_' + ed

                    geneID = ''

                    if spe1_prefix == 'sorghum':

                        print(col[8].split(';')[0])
                        if re.match('ID=.+_pg\d+\..+', col[8].split(';')[0]):
                            mt = re.match('ID=(.+_pg\d+)\..+', col[8].split(';')[0])
                            geneID = mt.group(1)
                        else:

                            mt = re.match('ID=(.+)', col[8].split(';')[0])
                            geneID = mt.group(1)

                    else:

                        if spe1_prefix == 'Uf':
                            if re.match('Name=(.+)', col[8].split(';')[1]):
                                mt = re.match('Name=(.+)', col[8].split(';')[1])
                                geneID = mt.group(1)

                        else:

                            if spe1_prefix == 'rice':

                                if re.match('Name=(.+)', col[8].split(';')[1]):
                                    mt = re.match('Name=(.+)', col[8].split(';')[1])
                                    geneID = mt.group(1)

                            else:
                                if re.match('ID=(.+)', col[8].split(';')[0]):
                                    mt = re.match('ID=(.+)', col[8].split(';')[0])
                                    geneID = mt.group(1)

                    if geneID != '':
                        #print(geneloc)
                        store_geneloc_geneID_dic[geneloc] = geneID

    store_acrloc_geneloc_dic = {}
    with open (ipt_spe1_acr_sumit_gene_closest_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[3]
            geneloc = col[7]
            if acrloc != '.' and geneloc != '.':
                store_acrloc_geneloc_dic[acrloc] = geneloc

    if spe1_prefix == 'rice':
        store_geneID_geneAnnot_dic = {}
        with open (ipt_spe1_gene_annot_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split('\t')
                geneID = col[0]
                if len(col) > 1:
                    annot = col[1]
                    store_geneID_geneAnnot_dic[geneID] = annot

    else:
        store_geneID_geneAnnot_dic = {}

    store_final_line_list = []
    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                acrloc = col[1]
                if acrloc in store_acrloc_geneloc_dic:
                    geneloc = store_acrloc_geneloc_dic[acrloc]
                    #print(geneloc)

                    geneID = store_geneloc_geneID_dic[geneloc]

                    if geneID in store_geneID_geneAnnot_dic:
                        annot = store_geneID_geneAnnot_dic[geneID]
                    else:
                        annot = 'none'
                else:
                    annot = 'none'
                    geneID = 'none'

                final_line = eachline + '\t' + geneID + '\t' + annot
                store_final_line_list.append(final_line)

            else:
                eachline = eachline + '\t' + 'geneID' + '\t' + 'Annot'
                store_final_line_list.append(eachline)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




##this will be copied to the s2_subfunctions.py
def subfunction_add_marker_gene_information (ipt_marker_gene_fl,ipt_final_add_gene_fl,opt_dir,spe1_prefix,spe2_prefix):


    store_marker_gene_celltype_dic = {}
    count = 0
    with open (ipt_marker_gene_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1
            if count != 1:
                geneID = col[0]
                celltype = col[2]
                tissue = col[1]
                organct = tissue + '__' + celltype
                store_marker_gene_celltype_dic[geneID] = organct


    store_final_line_list = []
    count = 0
    with open (ipt_final_add_gene_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            geneID = col[19]
            count += 1
            if count != 1:
                if geneID in store_marker_gene_celltype_dic:
                    organct = store_marker_gene_celltype_dic[geneID]
                else:
                    organct = 'none'

                final_line = eachline + '\t' + organct
                store_final_line_list.append(final_line)

            else:
                final_line = eachline + '\t' + 'Marker'
                store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_add_marker.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


##updating 120823
##we will add the gene expression data from the snRNA-seq rice
def subfunction_add_gene_exp_data_per_celltype (ipt_cpm_snRNAseq_fl,ipt_cluster_annot_fl,ipt_final_summary_with_geneID_fl,s2_s3_target_celltype_addCpm_str,opt_dir,
                                                spe1_prefix,spe2_prefix):

    ##/scratch/hy17471/rice_altlas_scATAC_seq_042021/05_integrate_snRNAseq_041322/01_Data_process_090922/04_annotate_by_scATACseq_data_042623/output_dir/step02_annotate_by_correlation_acc_exp_dir/opt_seedlingsnRNAseq_perM_genes_exp_clusters.txt
    ##/scratch/hy17471/rice_altlas_scATAC_seq_042021/05_integrate_snRNAseq_041322/01_Data_process_090922/04_annotate_by_scATACseq_data_042623/output_dir/step02_annotate_by_correlation_acc_exp_dir/opt_seedlingsnRNAseq_perM_genes_exp_clusters.txt

    store_cluster_celltype_dic = {}
    with open (ipt_cluster_annot_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            clusterID = col[0]
            annot = col[2]
            store_cluster_celltype_dic[clusterID] = annot

    store_gene_celltype_cpm_dic = {}
    store_order_clusterID_dic = {}
    count = 0
    with open (ipt_cpm_snRNAseq_fl,'r') as ipt:
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
                for i in range(1,len(col)):

                    celltype = store_order_clusterID_dic[i]
                    cpmval = col[i]

                    if celltype in s2_s3_target_celltype_addCpm_str.split(','):

                        store_celltype_cpmval_dic[celltype] = cpmval

                store_gene_celltype_cpm_dic[geneID] = store_celltype_cpmval_dic

    store_final_line_list = []
    count = 0
    with open (ipt_final_summary_with_geneID_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:

                geneID = col[19]

                if geneID in store_gene_celltype_cpm_dic:

                    store_celltype_cpmval_dic = store_gene_celltype_cpm_dic[geneID]
                    cpm_list = []
                    for eachcelltype in store_celltype_cpmval_dic:
                        cpm_list.append(store_celltype_cpmval_dic[eachcelltype])
                    cpm_line = '\t'.join(cpm_list)

                    final_line = eachline + '\t' + cpm_line
                    store_final_line_list.append(final_line)

                else:

                    none_list = []
                    for i in range(len(s2_s3_target_celltype_addCpm_str.split(','))):
                        none_list.append('none')

                    final_line = eachline + '\t' + '\t'.join(none_list)
                    store_final_line_list.append(final_line)

            else:

                celltype_list = []
                store_celltype_cpmval_dic = store_gene_celltype_cpm_dic[list(store_gene_celltype_cpm_dic.keys())[0]]
                for eachcelltype in store_celltype_cpmval_dic:
                    celltype_list.append(eachcelltype)

                final_line = eachline + '\t' + '\t'.join(celltype_list)
                store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_add_celltypeCpm.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



def subfunction_add_distance_to_celltypeCpm (ipt_summary_celltypeCpm_fl, opt_dir, ipt_gene_cate_fl,
                                             spe1_prefix,spe2_prefix):

    store_acr_gene_cate_dic = {}
    with open (ipt_gene_cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0]
            genecate = col[1]
            dist = col[2]

            store_acr_gene_cate_dic[acrloc] = genecate + '\t' + dist

    store_final_line_list = []
    count = 0
    with open (ipt_summary_celltypeCpm_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                acrloc = col[1]
                if acrloc in store_acr_gene_cate_dic:
                    final_line = eachline + '\t' + store_acr_gene_cate_dic[acrloc]
                else:
                    final_line = eachline + '\t' + 'none' + '\t' + 'none'

            else:
                final_line = eachline + '\t' + 'geneCate' + '\t' + 'Dist'

            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_add_celltypeCpm_addGeneCate.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




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




def multi_run(args):
    return subfunction_check_motif_enrichment_per_file(*args)


def subfunction_check_motif_enrichment_per_celltype (ipt_motif_fl,ipt_final_summary_fl,ipt_total_spe_acr_fl,ipt_opt_dir,
                                                     s2_s4_repeat_times,input_core_num):

    ##1. we will divide them into three category
    ##2. for each category, how many ACR per motif
    ##3. use the binomial test to check the enrichment

    store_shared_A_spe1_loc_dic = {}
    store_shared_I_spe1_celltype_loc_dic = {}
    store_not_shared_spe1_celltype_loc_dic = {}
    count = 0

    with open (ipt_final_summary_fl,'r') as ipt:
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

                ##for the shared-A
                if spe1_acrloc != 'none' and spe2_acrloc != 'none':
                    store_shared_A_spe1_loc_dic[spe1_acrloc] = 1

                ##for the shared-I
                if spe1_acrID != 'none' and spe2_acrloc == 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_shared_I_spe1_celltype_loc_dic:
                                store_shared_I_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_shared_I_spe1_celltype_loc_dic[eachspe1_celltype] = {}
                                store_shared_I_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                ##for the not shared
                if spe1_acrID == 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_not_shared_spe1_celltype_loc_dic:
                                store_not_shared_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_not_shared_spe1_celltype_loc_dic[eachspe1_celltype] = {}
                                store_not_shared_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

    ##Here we will save the multiple files into a dir
    store_all_target_acrloc_dir = ipt_opt_dir + '/store_all_target_acrloc_dir'
    if not os.path.exists(store_all_target_acrloc_dir):
        os.makedirs(store_all_target_acrloc_dir)

    two_category_list = ['shareI', 'notshare']
    for eachcate in two_category_list:

        target_loc_dic = {}
        if eachcate == 'shareI':
            target_loc_dic = store_shared_I_spe1_celltype_loc_dic
        if eachcate == 'notshare':
            target_loc_dic = store_not_shared_spe1_celltype_loc_dic

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
    pool.map(multi_run, run_list)


def subfunction_check_motif_enrichment_per_celltype_sharedA (ipt_motif_fl,ipt_final_summary_fl,ipt_total_spe_acr_fl,ipt_opt_dir,
                                                     s2_s4_repeat_times,input_core_num):

    ##1. we will divide them into three category
    ##2. for each category, how many ACR per motif
    ##3. use the binomial test to check the enrichment

    store_shared_A_spe1_loc_dic = {}
    store_shared_I_spe1_celltype_loc_dic = {}
    store_not_shared_spe1_celltype_loc_dic = {}
    store_shared_A_spe1_celltype_loc_dic = {}
    count = 0

    with open (ipt_final_summary_fl,'r') as ipt:
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

                ##for the shared-A
                if spe1_acrloc != 'none' and spe2_acrloc != 'none':
                    store_shared_A_spe1_loc_dic[spe1_acrloc] = 1

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_shared_A_spe1_celltype_loc_dic:
                                store_shared_A_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_shared_A_spe1_celltype_loc_dic[eachspe1_celltype] = {}
                                store_shared_A_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1


                ##for the shared-I
                if spe1_acrID != 'none' and spe2_acrloc == 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_shared_I_spe1_celltype_loc_dic:
                                store_shared_I_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_shared_I_spe1_celltype_loc_dic[eachspe1_celltype] = {}
                                store_shared_I_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

                ##for the not shared
                if spe1_acrID == 'none':

                    if spe1_celltype != 'broadly_accessible':
                        spe1_celltype_list = spe1_celltype.split(',')

                        for eachspe1_celltype in spe1_celltype_list:

                            if eachspe1_celltype in store_not_shared_spe1_celltype_loc_dic:
                                store_not_shared_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1
                            else:
                                store_not_shared_spe1_celltype_loc_dic[eachspe1_celltype] = {}
                                store_not_shared_spe1_celltype_loc_dic[eachspe1_celltype][spe1_acrloc] = 1

    ##Here we will save the multiple files into a dir
    store_all_target_acrloc_dir = ipt_opt_dir + '/store_all_target_acrloc_dir'
    if not os.path.exists(store_all_target_acrloc_dir):
        os.makedirs(store_all_target_acrloc_dir)

    two_category_list = ['shareA']
    for eachcate in two_category_list:

        target_loc_dic = {}
        if eachcate == 'shareA':
            target_loc_dic = store_shared_A_spe1_celltype_loc_dic

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
    pool.map(multi_run, run_list)



def subfunction_collect_all_enrich_results (ipt_spe_dir,opt_dir, spe1_prefix, spe2_prefix):


    all_temp_dir_list = glob.glob(ipt_spe_dir + '/store_parallele_analysis_dir/*')

    store_final_line_list = []
    for eachdir in all_temp_dir_list:

        final_enrich_fl = eachdir + '/opt_final_cate_motif_enrichment.txt'

        with open (final_enrich_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                final_line = eachline + '\t' + spe1_prefix + '_' + spe2_prefix
                store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_celltype_enrich.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


##updating 121023
def subfunction_identify_target_motif_binding_plotting (ipt_final_summary_all_spe_fl_dir,
                                                        ipt_motif_spe1_fimo_fl, ipt_target_motif_id_str,opt_dir,
                                                        spe1_prefix):

    ##we will first extract the target modif fimo location
    store_target_motif_ID_fimo_dir = opt_dir + '/store_target_motif_ID_fimo_dir'
    if not os.path.exists(store_target_motif_ID_fimo_dir):
        os.makedirs(store_target_motif_ID_fimo_dir)

    ipt_target_motif_id_str_list = ipt_target_motif_id_str.split(',')

    for eachmotifIDmotifnm in ipt_target_motif_id_str_list:

        mt = re.match('(.+):(.+)',eachmotifIDmotifnm)
        target_motifID = mt.group(1)
        target_motifnm = mt.group(2)

        store_final_line_list = []
        with open (ipt_motif_spe1_fimo_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                mt = re.match('(.+)_.+',col[3])
                motifID = mt.group(1)

                if motifID == target_motifID:
                    store_final_line_list.append(eachline)

        with open (store_target_motif_ID_fimo_dir + '/opt_' + target_motifnm + '.txt','w+') as opt:
            for eachline in store_final_line_list:
                opt.write(eachline + '\n')


    ##we will first check the syntenic regions across all species
    all_spe_final_summary_fl_list = glob.glob(ipt_final_summary_all_spe_fl_dir + '/*')
    store_spe_acr_dic = {}
    store_syntenic_reigon_notshared_acr_dic = {}
    store_notshared_acr_syntenic_reigon_dic = {}
    for eachspe_final_fl in all_spe_final_summary_fl_list:

        mt = re.match('.+/(.+)',eachspe_final_fl)
        flnm = mt.group(1)
        mt = re.match('(.+)_summary\.txt',flnm)
        spenm = mt.group(1)

        store_not_shared_ACR_dic = {}
        count = 0
        with open (eachspe_final_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                spe1acr_loc = col[1]
                count += 1
                if count != 1:

                    if col[0] == 'none':
                        store_not_shared_ACR_dic[spe1acr_loc] = 1
                        syntenic_region = col[11]

                        syntenic_region_list = syntenic_region.split(',')

                        for eachsyntenic_region in syntenic_region_list:

                            if eachsyntenic_region in store_syntenic_reigon_notshared_acr_dic:
                                store_syntenic_reigon_notshared_acr_dic[eachsyntenic_region][spe1acr_loc] = 1
                            else:
                                store_syntenic_reigon_notshared_acr_dic[eachsyntenic_region] = {}
                                store_syntenic_reigon_notshared_acr_dic[eachsyntenic_region][spe1acr_loc] = 1


                            if spe1acr_loc in store_notshared_acr_syntenic_reigon_dic:
                                store_notshared_acr_syntenic_reigon_dic[spe1acr_loc][eachsyntenic_region] = 1
                            else:
                                store_notshared_acr_syntenic_reigon_dic[spe1acr_loc] = {}
                                store_notshared_acr_syntenic_reigon_dic[spe1acr_loc][eachsyntenic_region] = 1

        store_spe_acr_dic[spenm] = store_not_shared_ACR_dic

    store_notshared_acr_spe_count_dic = {}
    for eachspe in store_spe_acr_dic:
        for eachacr in store_spe_acr_dic[eachspe]:
            if eachacr in store_notshared_acr_spe_count_dic:
                store_notshared_acr_spe_count_dic[eachacr] += 1
            else:
                store_notshared_acr_spe_count_dic[eachacr] = 1

    ##extract the acr shown in all species
    store_not_shared_ACR_inall_spe_dic = {}
    for eachacr in store_notshared_acr_spe_count_dic:

        if store_notshared_acr_spe_count_dic[eachacr] == 4:
            store_not_shared_ACR_inall_spe_dic[eachacr] = 1

    with open (opt_dir + '/temp_not_shared_acr_in_all_spe.txt','w+') as opt:
        for eachacr in store_not_shared_ACR_inall_spe_dic:
            acr_line = '\t'.join(eachacr.split('_'))
            opt.write(acr_line + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_not_shared_acr_in_all_spe.txt > ' + \
          opt_dir + '/temp_not_shared_acr_in_all_spe_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    all_motif_fl_list = glob.glob(store_target_motif_ID_fimo_dir + '/opt_*')
    for eachmotiffl in all_motif_fl_list:
        mt = re.match('.+/(.+)',eachmotiffl)
        flnm = mt.group(1)
        mt = re.match('opt_(.+)\.txt',flnm)
        motifnm = mt.group(1)

        cmd = 'bedtools intersect -wa -wb -a ' +  opt_dir + '/temp_not_shared_acr_in_all_spe_sorted.txt' + \
              ' -b ' + eachmotiffl + ' > ' + store_target_motif_ID_fimo_dir + '/temp_intersect_' + motifnm + '.txt'
        print(cmd)
        subprocess.call(cmd,shell=True)

    ##identify acr covering the number of target motifs
    store_motifnm_acr_motifloc_dic = {}
    all_intersect_motif_fl_list = glob.glob(store_target_motif_ID_fimo_dir + '/temp_intersect_*')
    for eachintersectmotif_fl in all_intersect_motif_fl_list:
        mt = re.match('.+/(.+)', eachintersectmotif_fl)
        flnm = mt.group(1)
        mt = re.match('temp_intersect_(.+)\.txt', flnm)
        motifnm = mt.group(1)

        store_acr_motifloc_dic = {}
        with open (eachintersectmotif_fl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()

                motifloc = col[3] + '_' + col[4] + '_' + col[5]
                acrloc = col[0] + '_' + col[1] + '_' + col[2]

                if acrloc in store_acr_motifloc_dic:
                    store_acr_motifloc_dic[acrloc][motifloc] = 1
                else:
                    store_acr_motifloc_dic[acrloc] = {}
                    store_acr_motifloc_dic[acrloc][motifloc] = 1

        store_motifnm_acr_motifloc_dic[motifnm] = store_acr_motifloc_dic

    ##we will build a summary to show all the acr and syntenic region and motif loc information
    store_final_line_list = []
    first_line = 'motif' + '\t' + 'acrloc' + '\t' + 'motifloc_str' + '\t' + 'motif_count' + '\t' + 'syntenic_region' + '\t' + 'acr_count_in_syntenic_region'
    store_final_line_list.append(first_line)
    for eachmotifnm in store_motifnm_acr_motifloc_dic:


        for eachacrloc in store_motifnm_acr_motifloc_dic[eachmotifnm]:

            motifloc_str = ','.join(list(store_motifnm_acr_motifloc_dic[eachmotifnm][eachacrloc].keys()))
            motifloc_count = len(list(store_motifnm_acr_motifloc_dic[eachmotifnm][eachacrloc].keys()))

            ##check the syntenic
            syntenic_region_dic = store_notshared_acr_syntenic_reigon_dic[eachacrloc]
            for eachsyntenic_region in syntenic_region_dic:

                acr_loc_within_syntenic_region_dic = store_syntenic_reigon_notshared_acr_dic[eachsyntenic_region]
                acr_loc_num = len(acr_loc_within_syntenic_region_dic.keys())

                final_line = eachmotifnm + '\t' + eachacrloc + '\t' + motifloc_str + '\t' + str(motifloc_count) + '\t'  + eachsyntenic_region + '\t' + str(acr_loc_num)
                store_final_line_list.append(final_line)




    with open (opt_dir + '/opt_final_' + spe1_prefix + '_target_motif_loc_in_not_shared_acr.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




#################
##udpating 122723
##we will add the two spe cpm other than only the rice for the variable ACR case
def subfunction_make_region_close_to_gene_cate (ipt_spe_gff, ipt_final_summary_fl,opt_dir,ipt_prefix):

    store_region_line_list = []
    count = 0
    with open(ipt_final_summary_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                spe2_region = col[4]

                if spe2_region != 'none':

                    spe2_region_list = spe2_region.split('_')
                    sumit = int((int(spe2_region_list[1]) + int(spe2_region_list[2])) / 2)

                    final_line = spe2_region_list[0] + '\t' + str(sumit) + '\t' + str(sumit + 1) + '\t' + spe2_region
                    # final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[0] + '_' + col[1] + '_' + col[2]
                    store_region_line_list.append(final_line)

    with open(opt_dir + '/temp_' + ipt_prefix + '_region_sumit_bed.txt', 'w+') as opt:
        for eachline in store_region_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_prefix + '_region_sumit_bed.txt > ' + \
          opt_dir + '/temp_' + ipt_prefix + '_region_sumit_bed_sorted.txt'
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
    cmd = 'bedtools closest -a ' + opt_dir + '/temp_' + ipt_prefix + '_region_sumit_bed_sorted.txt'+ \
          ' -b ' + opt_dir + '/temp_' + ipt_prefix + '_gene_bed_sorted.bed' + \
          ' -d > ' + opt_dir + '/temp_' + ipt_prefix + '_region_sumit_gene_closest.txt'
    subprocess.call(cmd, shell=True)

def subfunction_store_spe_gff_acrgeneloc_dic (ipt_spe_gff_fl,ipt_spe_acr_sumit_gene_closest_fl,spe_prefix):

    store_geneloc_geneID_dic = {}
    with open(ipt_spe_gff_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split()
                if col[2] == 'gene':

                    chrnm = col[0]
                    st = col[3]
                    ed = col[4]

                    geneloc = chrnm + '_' + st + '_' + ed

                    geneID = ''

                    if spe_prefix == 'sorghum':

                        print(col[8].split(';')[0])
                        if re.match('ID=.+_pg\d+\..+', col[8].split(';')[0]):
                            mt = re.match('ID=(.+_pg\d+)\..+', col[8].split(';')[0])
                            geneID = mt.group(1)
                        else:

                            mt = re.match('ID=(.+)', col[8].split(';')[0])
                            geneID = mt.group(1)

                    else:

                        if spe_prefix == 'Uf':
                            if re.match('Name=(.+)', col[8].split(';')[1]):
                                mt = re.match('Name=(.+)', col[8].split(';')[1])
                                geneID = mt.group(1)

                        else:

                            if spe_prefix == 'rice':

                                if re.match('Name=(.+)', col[8].split(';')[1]):
                                    mt = re.match('Name=(.+)', col[8].split(';')[1])
                                    geneID = mt.group(1)

                            else:
                                if re.match('ID=(.+)', col[8].split(';')[0]):
                                    mt = re.match('ID=(.+)', col[8].split(';')[0])
                                    geneID = mt.group(1)

                    if geneID != '':
                        store_geneloc_geneID_dic[geneloc] = geneID

    store_acrloc_geneloc_dic = {}
    with open (ipt_spe_acr_sumit_gene_closest_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[3]
            geneloc = col[7]
            store_acrloc_geneloc_dic[acrloc] = geneloc

    return (store_geneloc_geneID_dic,store_acrloc_geneloc_dic)


def subfunction_check_nearby_genes_two_spe (ipt_final_summary_fl,ipt_spe1_gff_fl,ipt_spe1_acr_sumit_gene_closest_fl,ipt_spe1_gene_annot_fl,
                                            ipt_spe2_gff_fl,ipt_spe2_region_gene_closest_fl,ipt_spe1_geneID_symbol_fl,
                                            opt_dir,spe1_prefix,spe2_prefix):


    store_geneloc_geneID_spe1_dic,store_acrloc_geneloc_spe1_dic = subfunction_store_spe_gff_acrgeneloc_dic (ipt_spe1_gff_fl,ipt_spe1_acr_sumit_gene_closest_fl,spe1_prefix)
    store_geneloc_geneID_spe2_dic,store_region_geneloc_spe2_dic = subfunction_store_spe_gff_acrgeneloc_dic(ipt_spe2_gff_fl, ipt_spe2_region_gene_closest_fl,spe2_prefix)

    store_geneID_symbol_spe1_dic = {}
    with open (ipt_spe1_geneID_symbol_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            geneID = col[0]
            symbol = col[1]
            store_geneID_symbol_spe1_dic[geneID] = symbol

    store_geneID_geneAnnot_spe1_dic = {}
    with open (ipt_spe1_gene_annot_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            geneID = col[0]
            if len(col) > 1:
                annot = col[1]
                store_geneID_geneAnnot_spe1_dic[geneID] = annot

    store_final_line_list = []
    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                acrloc_spe1 = col[1]

                if acrloc_spe1 in store_acrloc_geneloc_spe1_dic:
                    geneloc_spe1 = store_acrloc_geneloc_spe1_dic[acrloc_spe1]
                    geneID_spe1 = store_geneloc_geneID_spe1_dic[geneloc_spe1]

                    if geneID_spe1 in store_geneID_geneAnnot_spe1_dic:
                        annot_spe1 = store_geneID_geneAnnot_spe1_dic[geneID_spe1]
                    else:
                        annot_spe1 = 'none'

                    if geneID_spe1 in store_geneID_symbol_spe1_dic:
                        symbol_spe1 = store_geneID_symbol_spe1_dic[geneID_spe1]
                    else:
                        symbol_spe1 = 'none'

                else:
                    annot_spe1 = 'none'
                    geneID_spe1 = 'none'
                    symbol_spe1 = 'none'

                region_spe2 = col[4]
                if region_spe2 in store_region_geneloc_spe2_dic:
                    region_spe1 = store_region_geneloc_spe2_dic[region_spe2]
                    geneID_spe2 = store_geneloc_geneID_spe2_dic[region_spe1]

                else:
                    geneID_spe2 = 'none'

                final_line = eachline + '\t' + geneID_spe1 + '\t' + annot_spe1 + '\t' + symbol_spe1 + '\t' + geneID_spe2
                store_final_line_list.append(final_line)

            else:
                eachline = eachline + '\t' + spe1_prefix + '_geneID' + '\t' + spe1_prefix + '_Annot' + '\t' + spe1_prefix + '_Symbol' + '\t' + spe2_prefix + '_geneID'
                store_final_line_list.append(eachline)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_two_Spe.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




def subfunction_add_gene_exp_data_per_celltype_two_spe (ipt_cpm_snRNAseq_spe1_fl,ipt_cluster_annot_fl,ipt_cpm_snRNAseq_spe2_fl,ipt_final_summary_with_geneID_two_spe_fl,
                                                        s2_s3_target_celltype_addCpm_spe1_str,s2_s3_target_celltype_addCpm_spe2_str,opt_dir,spe1_prefix,spe2_prefix):

    ##/scratch/hy17471/rice_altlas_scATAC_seq_042021/05_integrate_snRNAseq_041322/01_Data_process_090922/04_annotate_by_scATACseq_data_042623/output_dir/step02_annotate_by_correlation_acc_exp_dir/opt_seedlingsnRNAseq_perM_genes_exp_clusters.txt
    ##/scratch/hy17471/rice_altlas_scATAC_seq_042021/05_integrate_snRNAseq_041322/01_Data_process_090922/04_annotate_by_scATACseq_data_042623/output_dir/step02_annotate_by_correlation_acc_exp_dir/opt_seedlingsnRNAseq_perM_genes_exp_clusters.txt

    store_cluster_celltype_dic = {}
    with open (ipt_cluster_annot_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            clusterID = col[0]
            annot = col[2]
            store_cluster_celltype_dic[clusterID] = annot

    store_gene_celltype_cpm_dic = {}
    store_order_clusterID_dic = {}
    count = 0
    with open (ipt_cpm_snRNAseq_spe1_fl,'r') as ipt:
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
                for i in range(1,len(col)):

                    celltype = store_order_clusterID_dic[i]
                    cpmval = col[i]

                    if celltype in s2_s3_target_celltype_addCpm_spe1_str.split(','):

                        store_celltype_cpmval_dic[celltype] = cpmval

                store_gene_celltype_cpm_dic[geneID] = store_celltype_cpmval_dic

    store_gene_celltype_cpm_spe2_dic = {}
    store_order_celltype_spe2_dic = {}
    count = 0
    with open (ipt_cpm_snRNAseq_spe2_fl,'r') as ipt:
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
                for i in range(1,len(col)):
                    celltype = store_order_celltype_spe2_dic[i]
                    cpmval = col[i]

                    if celltype in s2_s3_target_celltype_addCpm_spe2_str.split(','):
                        store_celltype_cpmval_dic[celltype] = cpmval



                store_gene_celltype_cpm_spe2_dic[geneID] = store_celltype_cpmval_dic

    #print(store_order_celltype_spe2_dic)
    #print(store_gene_celltype_cpm_spe2_dic)

    store_final_line_list = []
    count = 0
    with open (ipt_final_summary_with_geneID_two_spe_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1
            if count != 1:

                geneID_spe1 = col[19]

                if geneID_spe1 in store_gene_celltype_cpm_dic:

                    store_celltype_cpmval_dic = store_gene_celltype_cpm_dic[geneID_spe1]
                    cpm_list = []
                    for eachcelltype in store_celltype_cpmval_dic:
                        cpm_list.append(store_celltype_cpmval_dic[eachcelltype])
                    cpm_line_spe1 = '\t'.join(cpm_list)

                else:

                    none_list = []
                    for i in range(len(s2_s3_target_celltype_addCpm_spe1_str.split(','))):
                        none_list.append('none')

                    cpm_line_spe1 = '\t'.join(none_list)

                geneID_spe2 = col[22]
                if geneID_spe2 in store_gene_celltype_cpm_spe2_dic:

                    store_celltype_cpmval_dic = store_gene_celltype_cpm_spe2_dic[geneID_spe2]
                    cpm_list = []
                    for eachcelltype in store_celltype_cpmval_dic:
                        cpm_list.append(store_celltype_cpmval_dic[eachcelltype])
                    cpm_line_spe2 = '\t'.join(cpm_list)

                else:
                    none_list = []
                    for i in range(len(s2_s3_target_celltype_addCpm_spe2_str.split(','))):
                        none_list.append('none')

                    cpm_line_spe2 = '\t'.join(none_list)

                final_line = eachline + '\t' + cpm_line_spe1 + '\t' + cpm_line_spe2
                store_final_line_list.append(final_line)



            else:

                celltype_list = []
                store_celltype_cpmval_dic = store_gene_celltype_cpm_dic[list(store_gene_celltype_cpm_dic.keys())[0]]
                for eachcelltype in store_celltype_cpmval_dic:
                    celltype_list.append(eachcelltype)

                celltype_spe2_list = []
                store_celltype_cpmval_dic = store_gene_celltype_cpm_spe2_dic[list(store_gene_celltype_cpm_spe2_dic.keys())[0]]
                for eachcelltype in store_celltype_cpmval_dic:
                    celltype_spe2_list.append(eachcelltype)

                final_line = eachline + '\t' + '\t'.join(celltype_list) + '\t' + '\t'.join(celltype_spe2_list)
                store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_add_celltypeCpm_two_spe.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



#################
##udpating 122823
##we will add the two spe cpm other than only the rice for the species-specific case
def subfunction_add_gene_exp_data_per_celltype_two_spe_species_specificVer (ipt_final_summary_add_gene_ID_fl,ipt_cluster_annot_fl,ipt_cpm_snRNAseq_spe1_fl,
                                                                            ipt_cpm_snRNAseq_spe2_fl,ipt_ortho_fl,
                                                                            ipt_spe2_syntenic_region_gene_fl,opt_dir,
                                                                            s2_s3_target_celltype_addCpm_spe1_str,s2_s3_target_celltype_addCpm_spe2_str,
                                                                            spe1_prefix,spe2_prefix):

    store_spe1_spe2_geneID_dic = {}
    with open (ipt_ortho_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            mt = re.match('(.+)\.\d+',col[0])
            geneID = mt.group(1)
            store_spe1_spe2_geneID_dic[geneID] = col[1]

    store_cluster_celltype_dic = {}
    with open(ipt_cluster_annot_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            clusterID = col[0]
            annot = col[2]
            store_cluster_celltype_dic[clusterID] = annot

    store_gene_celltype_cpm_dic = {}
    store_order_clusterID_dic = {}
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

    store_gene_celltype_cpm_spe2_dic = {}
    store_order_celltype_spe2_dic = {}
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

    store_spe2_syntenic_region_gene_pair_list_dic = {}
    with open (ipt_spe2_syntenic_region_gene_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            syntenic_region = col[4]
            geneID = col[3]
            if syntenic_region in store_spe2_syntenic_region_gene_pair_list_dic:
                store_spe2_syntenic_region_gene_pair_list_dic[syntenic_region].append(geneID)
            else:
                store_spe2_syntenic_region_gene_pair_list_dic[syntenic_region] = []
                store_spe2_syntenic_region_gene_pair_list_dic[syntenic_region].append(geneID)

    store_final_line_list = []
    count = 0
    with open(ipt_final_summary_add_gene_ID_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')
            count += 1
            if count != 1:

                spe1_ACRID = col[0]
                spe1_celltype = col[8]
                spe1_ACR = col[1]

                if spe1_ACRID == 'none':

                    geneID_spe1 = col[19]

                    if geneID_spe1 in store_gene_celltype_cpm_dic:

                        store_celltype_cpmval_dic = store_gene_celltype_cpm_dic[geneID_spe1]
                        cpm_list = []
                        for eachcelltype in store_celltype_cpmval_dic:
                            cpm_list.append(store_celltype_cpmval_dic[eachcelltype])
                        cpm_line_spe1 = '\t'.join(cpm_list)

                    else:

                        none_list = []
                        for i in range(len(s2_s3_target_celltype_addCpm_spe1_str.split(','))):
                            none_list.append('none')

                        cpm_line_spe1 = '\t'.join(none_list)


                    syntenic_region = col[11]
                    if ',' not in syntenic_region:

                        spe2_gene_pair_list = store_spe2_syntenic_region_gene_pair_list_dic[syntenic_region]

                        for eachspe2gene in spe2_gene_pair_list:

                            geneID_spe2 = eachspe2gene
                            if geneID_spe2 in store_gene_celltype_cpm_spe2_dic:

                                store_celltype_cpmval_dic = store_gene_celltype_cpm_spe2_dic[geneID_spe2]
                                cpm_list = []
                                for eachcelltype in store_celltype_cpmval_dic:
                                    cpm_list.append(store_celltype_cpmval_dic[eachcelltype])
                                cpm_line_spe2 = '\t'.join(cpm_list)

                            else:
                                none_list = []
                                for i in range(len(s2_s3_target_celltype_addCpm_spe2_str.split(','))):
                                    none_list.append('none')

                                cpm_line_spe2 = '\t'.join(none_list)

                            if geneID_spe1 in store_spe1_spe2_geneID_dic:
                                orth_geneID_spe2 = store_spe1_spe2_geneID_dic[geneID_spe1]
                            else:
                                orth_geneID_spe2 = 'none'

                            if orth_geneID_spe2 == geneID_spe2:
                                final_orth = 'yes'
                            else:
                                final_orth = 'no'

                            final_line = spe1_ACR + '\t' + geneID_spe1 + '\t' + spe1_celltype + '\t' + cpm_line_spe1 + '\t' + syntenic_region + '\t' + geneID_spe2 + '\t' + cpm_line_spe2 + '\t' + final_orth
                            store_final_line_list.append(final_line)



            else:

                celltype_list = []
                store_celltype_cpmval_dic = store_gene_celltype_cpm_dic[list(store_gene_celltype_cpm_dic.keys())[0]]
                for eachcelltype in store_celltype_cpmval_dic:
                    celltype_list.append(eachcelltype)

                celltype_spe2_list = []
                store_celltype_cpmval_dic = store_gene_celltype_cpm_spe2_dic[list(store_gene_celltype_cpm_spe2_dic.keys())[0]]
                for eachcelltype in store_celltype_cpmval_dic:
                    celltype_spe2_list.append(eachcelltype)

                final_line = spe1_prefix + '_ACR' + '\t' + spe1_prefix + '_geneID' + '\t' + spe1_prefix + '_celltype' + '\t' + '\t'.join(celltype_list) + '\t' + 'syntenic_region' + '\t' + spe2_prefix + '_geneID' + '\t' + \
                             '\t'.join(celltype_spe2_list) + '\t' + 'IfRightOrth'
                store_final_line_list.append(final_line)


    with open(opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_add_celltypeCpm_two_spe_species_specificVer.txt',
            'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




















