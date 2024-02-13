#!/usr/bin/env python

##updating 112723 we will add the cross tissue analysis

import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats
from statistics import mean


def subfunction_check_all_maize_region_InOrNotH3K27me3 (ipt_target_organ_H3K27_peak_spe1_fl,ipt_target_organ_H3K27_peak_spe2_fl,
                                                        ipt_final_summary_blast_fl,
                                                        ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,
                                                        opt_dir,
                                                        s3_extend_H3K27me3_flank_bp,ipt_spe1_prefix,ipt_spe2_prefix):

    ##updating 111023
    ##add the rice gene cate
    store_acr_genecate_dist_spe1_dic = {}
    with open (ipt_spe1_acr_cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_acr_genecate_dist_spe1_dic[col[0]] = {'cate':col[1],'dist':col[2]}

    store_acr_genecate_dist_spe2_dic = {}
    with open (ipt_spe2_acr_cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_acr_genecate_dist_spe2_dic[col[0]] = {'cate':col[1],'dist':col[2]}


    store_maize_loc_dic = {}
    store_maize_acr_loc_dic = {}
    store_rice_acr_loc_dic = {}
    count = 0
    with open (ipt_final_summary_blast_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                maize_region = col[4]
                if maize_region != 'none':
                    maize_loc_line = '\t'.join(maize_region.split('_'))
                    store_maize_loc_dic[maize_loc_line] = 1

                    maize_acr = col[5]
                    if maize_acr != 'none':
                        maize_acr_loc_line = '\t'.join(maize_acr.split('_'))
                        store_maize_acr_loc_dic[maize_acr_loc_line] = 1

                rice_acr = col[1]
                if rice_acr != 'none':
                    rice_acr_loc_line = '\t'.join(rice_acr.split('_'))
                    store_rice_acr_loc_dic[rice_acr_loc_line] = 1



    with open (opt_dir + '/temp_' + ipt_spe2_prefix + '_blasted_region.txt','w+') as opt:
        for eachloc in store_maize_loc_dic:
            opt.write(eachloc + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_blasted_region.txt > ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_blasted_region_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    with open (opt_dir + '/temp_' + ipt_spe2_prefix + '_acr.txt','w+') as opt:
        for eachloc in store_maize_acr_loc_dic:
            opt.write(eachloc + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_acr.txt' + ' > ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    with open (opt_dir + '/temp_' + ipt_spe1_prefix + '_acr.txt','w+') as opt:
        for eachloc in store_rice_acr_loc_dic:
            opt.write(eachloc + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + ipt_spe1_prefix + '_acr.txt > ' + opt_dir + '/temp_' + ipt_spe1_prefix + '_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_final_line_list = []
    with open(ipt_target_organ_H3K27_peak_spe2_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            st = int(col[1]) - int(s3_extend_H3K27me3_flank_bp)
            if st < 0:
                new_st = 0
            else:
                new_st = st
            ed = int(col[2]) + int(s3_extend_H3K27me3_flank_bp)

            final_line = col[0] + '\t' + str(new_st) + '\t' + str(ed)
            store_final_line_list.append(final_line)

    with open(opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak.txt',
            'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak.txt' + ' > ' + opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##intersection
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_blasted_region_sorted.txt' + \
          ' -b ' + opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak_sorted.txt > ' + \
          opt_dir + '/temp_intersect_' + ipt_spe2_prefix + '_loc_' + 'H3K27me3_peak.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_maize_loc_blastedToH3K27me3_dic = {}
    with open (opt_dir + '/temp_intersect_' + ipt_spe2_prefix + '_loc_' + 'H3K27me3_peak.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            blastregion_loc = col[0] + '_' + col[1] + '_' + col[2]
            store_maize_loc_blastedToH3K27me3_dic[blastregion_loc] = 1

    ##intersect with acr
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + ipt_spe2_prefix + '_acr_sorted.txt' + \
          ' -b ' +  opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak_sorted.txt > ' + \
           opt_dir + '/temp_intersect_' + ipt_spe2_prefix + '_acr_' + 'H3K27me3_peak.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_maize_acrToH3K27me3_dic = {}
    with open (opt_dir + '/temp_intersect_' + ipt_spe2_prefix + '_acr_' + 'H3K27me3_peak.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            maize_acr = col[0] + '_' + col[1] + '_' + col[2]
            store_maize_acrToH3K27me3_dic[maize_acr] = 1


    ##intersect with acr in the rice
    store_final_line_list = []
    with open (ipt_target_organ_H3K27_peak_spe1_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            st = int(col[1]) - int(s3_extend_H3K27me3_flank_bp)
            if st < 0:
                new_st = 0
            else:
                new_st = st
            ed = int(col[2]) + int(s3_extend_H3K27me3_flank_bp)

            final_line = col[0] + '\t' + str(new_st) + '\t' + str(ed)
            store_final_line_list.append(final_line)

    with open(opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak_' + ipt_spe1_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + ipt_spe1_prefix + '_acr_sorted.txt' + \
          ' -b ' + opt_dir + '/temp_extended_' + s3_extend_H3K27me3_flank_bp + 'bp_H3K27m3_peak_' + ipt_spe1_prefix + '.txt > ' + \
          opt_dir + '/temp_intersect_' + ipt_spe1_prefix + '_acr_' + 'H3K27me3_peak.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_rice_acrToH3K27me3_dic = {}
    with open (opt_dir + '/temp_intersect_' + ipt_spe1_prefix + '_acr_' + 'H3K27me3_peak.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            rice_acr = col[0] + '_' + col[1] + '_' + col[2]
            store_rice_acrToH3K27me3_dic[rice_acr] = 1


    ##updating 110723
    ##Here we will build a new file updating if the ACR overlapping with the H3K27me3
    store_final_line_list = []
    count = 0
    with open(ipt_final_summary_blast_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                riceACR = col[1]
                maizeACR = col[5]
                maizeregion = col[4]
                if riceACR in store_rice_acrToH3K27me3_dic:
                    riceOverH3K27Cate = 'OverlapH3K27me3'
                else:
                    riceOverH3K27Cate = 'none'

                if maizeACR in store_maize_acrToH3K27me3_dic:
                    maizeOVerH3K27Cate = 'OverlapH3K27me3'
                else:
                    maizeOVerH3K27Cate = 'none'

                if riceACR in store_acr_genecate_dist_spe1_dic:
                    rice_genecate = store_acr_genecate_dist_spe1_dic[riceACR]['cate']
                    rice_genedist = store_acr_genecate_dist_spe1_dic[riceACR]['dist']
                else:
                    rice_genecate = 'none'
                    rice_genedist = 'none'

                if maizeACR in store_acr_genecate_dist_spe2_dic:
                    maize_genecate = store_acr_genecate_dist_spe2_dic[maizeACR]['cate']
                    maize_genedist = store_acr_genecate_dist_spe2_dic[maizeACR]['dist']

                else:
                    maize_genecate = 'none'
                    maize_genedist = 'none'

                if maizeregion in store_maize_loc_blastedToH3K27me3_dic:
                    maizeregionOverH3K27Cate = 'OverlapH3K27me3'
                else:
                    maizeregionOverH3K27Cate = 'none'

                final_line = eachline + '\t' + riceOverH3K27Cate + '\t' + maizeOVerH3K27Cate + '\t' + rice_genecate + '\t' + maize_genecate + '\t' + rice_genedist + '\t' + maize_genedist + '\t' + maizeregionOverH3K27Cate
                store_final_line_list.append(final_line)
            else:
                final_line = eachline + '\t' + ipt_spe1_prefix + 'OverH3K27Cate' + '\t' + ipt_spe2_prefix + 'OVerH3K27Cate' + '\t' + ipt_spe1_prefix + 'GeneCate' + '\t' + ipt_spe2_prefix + 'GeneCate' + \
                            '\t' + ipt_spe1_prefix + 'GeneDist' + '\t' + ipt_spe2_prefix + 'GeneDist' + '\t' + ipt_spe2_prefix + 'RegionOverH3K27Cate'
                store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_final_summary_file_add_' + ipt_spe1_prefix + '_' + ipt_spe2_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    return (store_maize_loc_blastedToH3K27me3_dic,store_maize_acrToH3K27me3_dic,store_rice_acrToH3K27me3_dic)


##updating 010324
##we will build a new function to add the species specific cases so there is no reigons overlapped with the H3K27me3
##we found it is okay if nothing has been found in the maize





def subfunction_store_different_cate_riceACRs (maizeACRloc,maize_H3K27me3Cate,maize_BroOrRes,store_maize_acrToH3K27me3_dic,
                                               riceACRloc,maizeblastedRg,store_maize_loc_blastedToH3K27me3_dic,
                                               store_broad_to_broad_riceACRnotCoverH3K27me3_dic,store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic,
                                               store_broad_to_noACR_notCoverH3K27me3_dic):



    ##There are something wrong for some

    ##for the broad to broad we will consider two cases under the H3K27me3 and not for the maize
    if maizeACRloc != 'none':

        if maize_H3K27me3Cate == 'BroadInflankH3K27me3peak':
            cate = 'H3K27me3'

            if cate in store_broad_to_broad_riceACRnotCoverH3K27me3_dic:
                store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1
            else:
                store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate] = {}
                store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1

        else:

            if maize_H3K27me3Cate != 'RestrictedInflankH3K27me3peak':

                ##make sure maizeACRloc not cover H3K27me3
                if maizeACRloc not in store_maize_acrToH3K27me3_dic:

                    if maize_BroOrRes == 'broad':
                        cate = 'H3K27me3Absent'

                        if cate in store_broad_to_broad_riceACRnotCoverH3K27me3_dic:
                            store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1
                        else:
                            store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate] = {}
                            store_broad_to_broad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1

                    else:

                        cate = 'H3K27me3Absent'
                        if cate in store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic:
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1
                        else:
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate] = {}
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1

                ##another cate to include the others underlying the H3K27me3
                else:

                    if maize_BroOrRes == 'others' or maize_BroOrRes == 'restricted':
                        cate = 'H3K27me3'
                        if cate in store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic:
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1
                        else:
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate] = {}
                            store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1

                    ##Here the else the broad we already consider the above condition


            else:
                ##Here is only for the restriced in H3K27me3peak we will also check the others in the H3K27me3
                cate = 'H3K27me3'
                if cate in store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic:
                    store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1
                else:
                    store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate] = {}
                    store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic[cate][riceACRloc] = 1


    ## if there is no acr
    else:
        if maizeblastedRg != 'none':

            ##for the H3K27me3
            if maizeblastedRg in store_maize_loc_blastedToH3K27me3_dic:
                cate = 'H3K27me3'
                if cate in store_broad_to_noACR_notCoverH3K27me3_dic:
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate][riceACRloc] = 1
                else:
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate] = {}
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate][riceACRloc] = 1
            else:
                cate = 'H3K27me3Absent'
                if cate in store_broad_to_noACR_notCoverH3K27me3_dic:
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate][riceACRloc] = 1
                else:
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate] = {}
                    store_broad_to_noACR_notCoverH3K27me3_dic[cate][riceACRloc] = 1



def subfunction_summarize_overview_for_H3K27me3 (ipt_final_summary_blast_fl,ipt_target_organ_H3K27_peak_spe1_fl,
                                                 ipt_target_organ_H3K27_peak_spe2_fl,
                                                 ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,
                                                 spe1_prefix,spe2_prefix,opt_dir,
                                                 s3_extend_H3K27me3_flank_bp):

    ##we need to first store blasted loc whether they are overlapped with the H3K27me3
    store_maize_loc_blastedToH3K27me3_dic, store_maize_acrToH3K27me3_dic,\
    store_rice_acrToH3K27me3_dic = subfunction_check_all_maize_region_InOrNotH3K27me3(ipt_target_organ_H3K27_peak_spe1_fl,
                                                                                      ipt_target_organ_H3K27_peak_spe2_fl,
                                                                                      ipt_final_summary_blast_fl,
                                                                                      ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,
                                                                                      opt_dir,
                                                                                      s3_extend_H3K27me3_flank_bp, spe1_prefix,spe2_prefix)

    ##Looks like the following may not work as we will use the output from the subfunction_check_all_maize_region_InOrNotH3K27me3
    ##to do the plot in the R

    #################
    ##category type 1

    ##H3K27me3 not shared
    count = 0
    store_H3K27me3_absent_riceACR_dic = {}
    store_H3K27me3_sharedbroad_acc_riceACR_dic = {}
    store_H3K27me3_sharedrestricted_acc_riceACR_dic = {}
    store_H3K27me3_sharedinacc_riceACR_dic = {}
    with open(ipt_final_summary_blast_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                riceACRloc = col[1]
                maizeACRloc = col[5]
                maizeblastedRg = col[4]
                rice_H3K27me3Cate = col[2]
                maize_H3K27me3Cate = col[6]

                if rice_H3K27me3Cate == 'BroadInflankH3K27me3peak':

                    if maizeACRloc != 'none':
                        ##indicate maize ACR not overlap with the H3K27me3
                        if maizeACRloc not in store_maize_acrToH3K27me3_dic:
                            store_H3K27me3_absent_riceACR_dic[riceACRloc] = 1

                        else:
                            if maize_H3K27me3Cate == 'BroadInflankH3K27me3peak':
                                store_H3K27me3_sharedbroad_acc_riceACR_dic[riceACRloc] = 1

                            if maize_H3K27me3Cate == 'RestrictedInflankH3K27me3peak':
                                store_H3K27me3_sharedrestricted_acc_riceACR_dic[riceACRloc] = 1

                    else:
                        if maizeblastedRg != 'none':
                            if maizeblastedRg not in store_maize_loc_blastedToH3K27me3_dic:
                                store_H3K27me3_absent_riceACR_dic[riceACRloc] = 1
                            else:
                                store_H3K27me3_sharedinacc_riceACR_dic[riceACRloc] = 1

    ##report the results
    store_final_line_list = []
    final_line = 'H3K27me3_absent' + '\t' + str(len(list(store_H3K27me3_absent_riceACR_dic.keys())))
    store_final_line_list.append(final_line)
    final_line = 'H3K27me3_shared_broad_accessible' + '\t' + str(len(list(store_H3K27me3_sharedbroad_acc_riceACR_dic.keys())))
    store_final_line_list.append(final_line)
    final_line = 'H3K27me3_shared_restricted_accessible' + '\t' + str(len(list(store_H3K27me3_sharedrestricted_acc_riceACR_dic.keys())))
    store_final_line_list.append(final_line)
    final_line = 'H3K27me3_shared_inaccessible' + '\t' + str(len(list(store_H3K27me3_sharedinacc_riceACR_dic.keys())))
    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_summary_H3K27me3_broad_' + spe1_prefix + 'ACR_to_' + spe2_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    #################
    ##category type 2
    ##This is for the covering H3K27me3
    store_broad_to_broad_riceACR_dic = {}
    store_broad_to_notbroad_riceACR_dic = {}
    store_broad_to_noACR_dic = {}

    ##updating 110623
    ##Here we will build a background to use the fisher to check the correspondings were enriched for the ACR under the H3K27me3
    store_broad_to_broad_riceACRnotCoverH3K27me3_dic = {}
    store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic = {}
    store_broad_to_noACR_notCoverH3K27me3_dic = {}

    store_notbroad_to_broad_riceACR_dic = {}
    store_notbroad_to_notbroad_riceACR_dic = {}
    store_notbroad_to_noACR_dic = {}

    with open(ipt_final_summary_blast_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                riceACRloc = col[1]
                rice_BroOrRes = col[3]
                maizeACRloc = col[5]
                maizeblastedRg = col[4]
                rice_H3K27me3Cate = col[2]
                maize_H3K27me3Cate = col[6]
                maize_BroOrRes = col[7]

                if rice_H3K27me3Cate == 'BroadInflankH3K27me3peak':

                    subfunction_store_different_cate_riceACRs(maizeACRloc, maize_H3K27me3Cate, maize_BroOrRes,
                                                              store_maize_acrToH3K27me3_dic,
                                                              riceACRloc, maizeblastedRg,
                                                              store_maize_loc_blastedToH3K27me3_dic,
                                                              store_broad_to_broad_riceACR_dic,
                                                              store_broad_to_notbroad_riceACR_dic,
                                                              store_broad_to_noACR_dic)


                #####################################################
                ##Here we will store the rice ACR not in the H3K27me3
                else:
                    if riceACRloc not in store_rice_acrToH3K27me3_dic:

                        if rice_BroOrRes == 'broad':

                            subfunction_store_different_cate_riceACRs(maizeACRloc, maize_H3K27me3Cate, maize_BroOrRes,
                                                                      store_maize_acrToH3K27me3_dic,
                                                                      riceACRloc, maizeblastedRg,
                                                                      store_maize_loc_blastedToH3K27me3_dic,
                                                                      store_broad_to_broad_riceACRnotCoverH3K27me3_dic,
                                                                      store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic,
                                                                      store_broad_to_noACR_notCoverH3K27me3_dic)


                    else:


                        subfunction_store_different_cate_riceACRs(maizeACRloc, maize_H3K27me3Cate, maize_BroOrRes,
                                                                  store_maize_acrToH3K27me3_dic,
                                                                  riceACRloc, maizeblastedRg,
                                                                  store_maize_loc_blastedToH3K27me3_dic,
                                                                  store_notbroad_to_broad_riceACR_dic,
                                                                  store_notbroad_to_notbroad_riceACR_dic,
                                                                  store_notbroad_to_noACR_dic)


    ##now we will check the proportion for all cate
    store_final_line_list = []

    total_broad_to_broad_riceACR_dic = {}
    for eachcate in store_broad_to_broad_riceACR_dic:
        for eachacr in store_broad_to_broad_riceACR_dic[eachcate]:
            total_broad_to_broad_riceACR_dic[eachacr] = 1
    total_broad_to_broad_riceACR_num = len(list(total_broad_to_broad_riceACR_dic.keys()))
    for eachcate in store_broad_to_broad_riceACR_dic:
        cate_acr_dic = store_broad_to_broad_riceACR_dic[eachcate]
        cate_acr_num = len(list(cate_acr_dic.keys()))
        final_line = 'BroadToBroad' + '\t' + eachcate + '\t' + str(cate_acr_num) + '\t' + str(total_broad_to_broad_riceACR_num)
        store_final_line_list.append(final_line)

    total_broad_to_notbroad_riceACR_dic = {}
    for eachcate in store_broad_to_notbroad_riceACR_dic:
        for eachacr in store_broad_to_notbroad_riceACR_dic[eachcate]:
            total_broad_to_notbroad_riceACR_dic[eachacr] = 1
    total_broad_to_notbroad_riceACR_num = len(list(total_broad_to_notbroad_riceACR_dic.keys()))
    for eachcate in store_broad_to_notbroad_riceACR_dic:
        cate_acr_dic = store_broad_to_notbroad_riceACR_dic[eachcate]
        cate_acr_num = len(list(cate_acr_dic.keys()))
        final_line = 'BroadToNotBroad' + '\t' + eachcate + '\t' + str(cate_acr_num) + '\t' + str(total_broad_to_notbroad_riceACR_num)
        store_final_line_list.append(final_line)

    total_broad_to_noACR_riceACR_dic = {}
    for eachcate in store_broad_to_noACR_dic:
        for eachacr in store_broad_to_noACR_dic[eachcate]:
            total_broad_to_noACR_riceACR_dic[eachacr] = 1
    total_broad_to_noACR_riceACR_num = len(list(total_broad_to_noACR_riceACR_dic.keys()))
    for eachcate in store_broad_to_noACR_dic:
        cate_acr_dic = store_broad_to_noACR_dic[eachcate]
        cate_acr_num = len(list(cate_acr_dic.keys()))
        final_line = 'BroadToNoACR' + '\t' + eachcate + '\t' + str(cate_acr_num) + '\t' + str(total_broad_to_noACR_riceACR_num)
        store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_summary_H3K27me3_broad_' + spe1_prefix + 'ACR_to_' + spe2_prefix + '_cateType_second_Version.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##updating 110623
    ##conduct the enrichment test
    ##See if ACR under H3K27me3 in rice would like to be enriched on the ACR under the H3K27me3
    ##for the broad to broad
    store_final_line_list = []

    BtoB_H3K27rice_H3K27maize = len(list(store_broad_to_broad_riceACR_dic['H3K27me3'].keys()))
    BtoB_H3K27rice_absentmaize = len(list(store_broad_to_broad_riceACR_dic['H3K27me3Absent'].keys()))

    BtoNB_H3K27rice_H3K27maize = len(list(store_broad_to_notbroad_riceACR_dic['H3K27me3'].keys()))
    BtoNB_H3K27rice_absentmaize = len(list(store_broad_to_notbroad_riceACR_dic['H3K27me3Absent'].keys()))

    BtoNoACR_H3K27rice_H3K27maize = len(list(store_broad_to_noACR_dic['H3K27me3'].keys()))
    BtoNoACR_H3K27rice_absentmaize = len(list(store_broad_to_noACR_dic['H3K27me3Absent'].keys()))

    ##for the ACR in rice not overlapped with the H3K27me3
    BtoB_notH3K27rice_H3K27maize = len(list(store_broad_to_broad_riceACRnotCoverH3K27me3_dic['H3K27me3'].keys()))
    BtoB_notH3K27rice_absentmaize = len(list(store_broad_to_broad_riceACRnotCoverH3K27me3_dic['H3K27me3Absent'].keys()))

    BtoNB_notH3K27rice_H3K27maize = len(list(store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic['H3K27me3'].keys()))
    BtoNB_notH3K27rice_absentmaize = len(list(store_broad_to_notbroad_riceACRnotCoverH3K27me3_dic['H3K27me3Absent'].keys()))

    BtoNoACR_notH3K27rice_H3K27maize = len(list(store_broad_to_noACR_notCoverH3K27me3_dic['H3K27me3'].keys()))
    BtoNoACR_notH3K27rice_absentmaize = len(list(store_broad_to_noACR_notCoverH3K27me3_dic['H3K27me3Absent'].keys()))

    ##for the notbroad rice ACR under H3K27me3 to all the other cases
    NBtoB_H3K27rice_H3K27maize = len(list(store_notbroad_to_broad_riceACR_dic['H3K27me3'].keys()))
    NBtoB_H3K27rice_absentmaize = len(list(store_notbroad_to_broad_riceACR_dic['H3K27me3Absent'].keys()))

    NBtoNB_H3K27rice_H3K27maize = len(list(store_notbroad_to_notbroad_riceACR_dic['H3K27me3'].keys()))
    NBtoNB_H3K27rice_absentmaize = len(list(store_notbroad_to_notbroad_riceACR_dic['H3K27me3Absent'].keys()))

    NBtoNoACR_H3K27rice_H3K27maize = len(list(store_notbroad_to_noACR_dic['H3K27me3'].keys()))
    NBtoNoACR_H3K27rice_absentmaize = len(list(store_notbroad_to_noACR_dic['H3K27me3Absent'].keys()))


    ##use the fisher to check the enrichment
    oddsratio, pvalue = stats.fisher_exact([[BtoB_H3K27rice_H3K27maize, BtoB_H3K27rice_absentmaize],
                                            [BtoB_notH3K27rice_H3K27maize, BtoB_notH3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoB' + '\t' + str(pvalue) + '\t' + str(BtoB_H3K27rice_H3K27maize) + '\t' + str(BtoB_H3K27rice_absentmaize) + '\t' + str(BtoB_notH3K27rice_H3K27maize) + '\t' + str(BtoB_notH3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    oddsratio, pvalue = stats.fisher_exact([[BtoNB_H3K27rice_H3K27maize, BtoNB_H3K27rice_absentmaize],
                                            [BtoNB_notH3K27rice_H3K27maize, BtoNB_notH3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoNB' + '\t' + str(pvalue) + '\t' + str(BtoNB_H3K27rice_H3K27maize) + '\t' + str(
        BtoNB_H3K27rice_absentmaize) + '\t' + str(BtoNB_notH3K27rice_H3K27maize) + '\t' + str(
        BtoNB_notH3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    oddsratio, pvalue = stats.fisher_exact([[BtoNoACR_H3K27rice_H3K27maize, BtoNoACR_H3K27rice_absentmaize],
                                            [BtoNoACR_notH3K27rice_H3K27maize, BtoNoACR_notH3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoNoACR' + '\t' + str(pvalue) + '\t' + str(BtoNoACR_H3K27rice_H3K27maize) + '\t' + str(
        BtoNoACR_H3K27rice_absentmaize) + '\t' + str(BtoNoACR_notH3K27rice_H3K27maize) + '\t' + str(
        BtoNoACR_notH3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_summary_H3K27me3_broad_' + spe1_prefix + 'ACR_to_' + spe2_prefix + '_cateType_second_Version_fisherTest_CompareTonotH3K27rice.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    store_final_line_list = []
    oddsratio, pvalue = stats.fisher_exact([[BtoB_H3K27rice_H3K27maize, BtoB_H3K27rice_absentmaize],
                                            [NBtoB_H3K27rice_H3K27maize, NBtoB_H3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoB' + '\t' + str(pvalue) + '\t' + str(BtoB_H3K27rice_H3K27maize) + '\t' + str(BtoB_H3K27rice_absentmaize) + '\t' + str(NBtoB_H3K27rice_H3K27maize) + '\t' + str(NBtoB_H3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    oddsratio, pvalue = stats.fisher_exact([[BtoNB_H3K27rice_H3K27maize, BtoNB_H3K27rice_absentmaize],
                                            [NBtoNB_H3K27rice_H3K27maize, NBtoNB_H3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoNB' + '\t' + str(pvalue) + '\t' + str(BtoNB_H3K27rice_H3K27maize) + '\t' + str(
        BtoNB_H3K27rice_absentmaize) + '\t' + str(NBtoNB_H3K27rice_H3K27maize) + '\t' + str(
        NBtoNB_H3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    oddsratio, pvalue = stats.fisher_exact([[BtoNoACR_H3K27rice_H3K27maize, BtoNoACR_H3K27rice_absentmaize],
                                            [NBtoNoACR_H3K27rice_H3K27maize, NBtoNoACR_H3K27rice_absentmaize]],
                                           alternative='greater')
    final_line = 'BtoNoACR' + '\t' + str(pvalue) + '\t' + str(BtoNoACR_H3K27rice_H3K27maize) + '\t' + str(
        BtoNoACR_H3K27rice_absentmaize) + '\t' + str(NBtoNoACR_H3K27rice_H3K27maize) + '\t' + str(
        NBtoNoACR_H3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_summary_H3K27me3_broad_' + spe1_prefix + 'ACR_to_' + spe2_prefix + '_cateType_second_Version_fisherTest_CompareToNotBroadH3K27me3rice.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##check if H3K27me3 covered BtoB
    store_final_line_list = []
    oddsratio, pvalue = stats.fisher_exact([[BtoB_H3K27rice_H3K27maize, BtoNB_H3K27rice_H3K27maize],
                                            [BtoB_notH3K27rice_absentmaize, BtoNB_notH3K27rice_absentmaize]],
                                           alternative='less')

    final_line = 'H3K27compareNotH3K27me3' + '\t' + str(pvalue) + '\t' + str(BtoB_H3K27rice_H3K27maize) + '\t' + str(BtoNB_H3K27rice_H3K27maize) + '\t' + str(BtoB_notH3K27rice_absentmaize) + '\t' + \
        str(BtoNB_notH3K27rice_absentmaize)
    store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_summary_H3K27me3_broad_' + spe1_prefix + 'ACR_to_' + spe2_prefix + '_cateType_second_Version_fisherTest_H3K27compareNotH3K27me3.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##use the fisher to check the enrichment of broad compare to other and restricted ones



##updating 112723
########
##step03
##want to check if the conserved H3K27me3 peaks were also conserved in other organs
def subfunction_check_conserved_H3K27me3_in_other_organs (ipt_final_summary_blast_fl,ipt_other_organ_fl,opt_dir,
                                                          spe1_prefix,spe2_prefix,organ_prefix):


    ##we will extract the target column
    store_final_line_list = []
    count = 0
    with open (ipt_final_summary_blast_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                ACR_line = '\t'.join(col[1].split('_'))

                target_line = col[2] + '\t' + col[3] + '\t' + col[4] + '\t' + col[5] + '\t' + col[6] + '\t' + col[7]
                final_line = ACR_line + '\t' + target_line
                store_final_line_list.append(final_line)

    with open (opt_dir + '/temp_final_summary_file_concise.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_final_summary_file_concise.txt > ' + \
          opt_dir + '/temp_final_summary_file_concise_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##we will also extract the essential information in the second part
    store_final_line_list = []
    count = 0
    with open (ipt_other_organ_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[9] + '\t' + col[10]
                store_final_line_list.append(final_line)

    with open (opt_dir + '/temp_final_' + organ_prefix + '_H3K27me3_fl.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_final_' + organ_prefix + '_H3K27me3_fl.txt > ' + \
          opt_dir + '/temp_final_other_organ_H3K27me3_fl_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##check the intersection
    cmd = 'bedtools intersect -wa -wb -a ' +  opt_dir + '/temp_final_summary_file_concise_sorted.txt' + \
          ' -b ' + opt_dir + '/temp_final_other_organ_H3K27me3_fl_sorted.txt > ' + opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_' + organ_prefix + '.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)


########
##step04
##updating 121223
##we will check if there are motif difference for the different cate
def subfunction_check_motif_freq_in_diff_cate (ipt_final_summary_fl, ipt_motif_fimo_fl,opt_dir, spe1_prefix,spe2_prefix):

    ##select the target ACRs
    store_shared_A_acr_dic = {}
    store_shared_I_acr_dic = {}
    store_not_shared_acr_dic = {}
    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                riceACRID = col[0]
                riceACRloc = col[1]
                maizeACRloc = col[5]
                riceInorH3K27me3 = col[2]

                ##for the shared A ACRs
                if riceACRID != 'none' and maizeACRloc != 'none':
                    if riceInorH3K27me3 == 'BroadInflankH3K27me3peak':
                        store_shared_A_acr_dic[riceACRloc] = 1

                ##for the shared I ACRs
                if riceACRID != 'none' and maizeACRloc == 'none':
                    if riceInorH3K27me3 == 'BroadInflankH3K27me3peak':
                        store_shared_I_acr_dic[riceACRloc] = 1

                if riceACRID == 'none':
                    if riceInorH3K27me3 == 'BroadInflankH3K27me3peak':
                        store_not_shared_acr_dic[riceACRloc] = 1


    with open (opt_dir + '/temp_all_cate_acr_loc_underH3K27.txt','w+') as opt:
        for eachacr in store_shared_A_acr_dic:
            acrline = '\t'.join(eachacr.split('_'))
            final_line = acrline + '\t' + 'sharedA'
            opt.write(final_line + '\n')
        for eachacr in store_shared_I_acr_dic:
            acrline = '\t'.join(eachacr.split('_'))
            final_line = acrline + '\t' + 'sharedI'
            opt.write(final_line + '\n')
        for eachacr in store_not_shared_acr_dic:
            acrline = '\t'.join(eachacr.split('_'))
            final_line = acrline + '\t' + 'not_shared'
            opt.write(final_line + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_all_cate_acr_loc_underH3K27.txt > ' + \
          opt_dir + '/temp_all_cate_acr_loc_underH3K27_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_all_cate_acr_loc_underH3K27_sorted.txt' + \
          ' -b ' + ipt_motif_fimo_fl + ' > ' + \
          opt_dir + '/temp_intersect.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    ##check the proportion of ACRs covering each motif
    store_motif_cate_acr_dic = {}
    with open (opt_dir + '/temp_intersect.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            mt = re.match('.+-(.+)_.+',col[7])
            motifnm = mt.group(1)

            cate = col[3]
            acrloc = col[0] + '_' + col[1] + '_' + col[2]

            if motifnm in store_motif_cate_acr_dic:
                if cate in store_motif_cate_acr_dic[motifnm]:
                    store_motif_cate_acr_dic[motifnm][cate][acrloc] = 1
                else:
                    store_motif_cate_acr_dic[motifnm][cate] = {}
                    store_motif_cate_acr_dic[motifnm][cate][acrloc] = 1
            else:
                store_motif_cate_acr_dic[motifnm] = {}
                store_motif_cate_acr_dic[motifnm][cate] = {}
                store_motif_cate_acr_dic[motifnm][cate][acrloc] = 1

    store_final_line_list = []
    for eachmotifnm in store_motif_cate_acr_dic:
        for eachcate in store_motif_cate_acr_dic[eachmotifnm]:

            acrloc_cover_motif_dic = store_motif_cate_acr_dic[eachmotifnm][eachcate]

            total_acrloc_dic = {}
            if eachcate == 'sharedA':
                total_acrloc_dic = store_shared_A_acr_dic
            else:
                if eachcate == 'sharedI':
                    total_acrloc_dic = store_shared_I_acr_dic
                else:
                    if eachcate == 'not_shared':
                        total_acrloc_dic = store_not_shared_acr_dic
                    else:
                        print(eachcate)

            if total_acrloc_dic != {}:

                total_acrloc_num = len(list(total_acrloc_dic.keys()))

                acrloc_cover_motif_num = len(list(acrloc_cover_motif_dic.keys()))

                final_line = eachmotifnm + '\t' + eachcate + '\t' + str(acrloc_cover_motif_num/total_acrloc_num) + '\t' + str(acrloc_cover_motif_num) + '\t' + str(total_acrloc_num)
                store_final_line_list.append(final_line)

            else:
                print('there is somethingf wrong with the eachcate')

    with open (opt_dir + '/opt_final_' + spe1_prefix + '_' + spe2_prefix + '_motif_cate_acr_num.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')




##updating 122023
##we will check the overlapping with the CNS
def subfunction_overlap_with_CNS (ipt_spe1_final_summary_fl,ipt_spe1_cns_gff_fl,spe1_prefix,spe2_prefix,opt_dir):


    ##build the acr file
    store_acr_line_dic = {}
    count = 0
    with open (ipt_spe1_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split()
                acr_line = '\t'.join(col[1].split('_'))
                store_acr_line_dic[acr_line] = 1

    with open (opt_dir + '/temp_acr_fl.txt','w+') as opt:
        for eachline in store_acr_line_dic:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_acr_fl.txt > ' + \
          opt_dir + '/temp_acr_fl_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

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
        with open(ipt_spe1_cns_gff_fl, 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                if not eachline.startswith('#'):
                    col = eachline.strip().split()
                    cns_line = col[0] + '\t' + col[3] + '\t' + col[4]
                    store_final_line_list.append(cns_line)

    else:
        store_final_line_list = []
        with open(ipt_spe1_cns_gff_fl, 'r') as ipt:
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

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_acr_fl.txt' + \
          ' -b ' + opt_dir + '/temp_cns_sorted.txt > ' + \
          opt_dir + '/temp_intersect_acr_cns.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_acr_overlapcns_dic = {}
    with open (opt_dir + '/temp_intersect_acr_cns.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[0] + '_' + col[1] + '_' + col[2]
            store_acr_overlapcns_dic[acrloc] = 1

    store_final_line_list = []
    count = 0
    with open (ipt_spe1_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                spe1_acr = col[1]
                if spe1_acr in store_acr_overlapcns_dic:
                    overlapCNS = 'yes'
                else:
                    overlapCNS = 'no'

                final_line = eachline + '\t' + overlapCNS
            else:
                final_line = eachline + '\t' + 'OverlapCNS'

            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_final_summary_file_add_CNS_' + spe1_prefix + '_' + spe2_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


def subfunction_overlap_with_atlas_acr (ipt_rice_atlas_fl,input_rice_atlas_acr_coverage_fl,ipt_final_summary_fl,
                                        ipt_rice_all_organ_H3K27me3_dir,
                                        s2_s2_organct_coverage_cutoff,opt_dir,spe1_prefix,spe2_prefix):


    ##save the cell type num for the spe2
    ##we will check the H3K27m3 for each of organ one by one
    store_all_rice_atlas_peak_dic = {}
    with open(ipt_rice_atlas_fl, 'r') as ipt:
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

    with open(opt_dir + '/temp_' + spe1_prefix + '_acr.txt', 'w+') as opt:
        for eachline in store_spe1_acr_dic:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe1_prefix + '_acr.txt > ' + \
          opt_dir + '/temp_' + spe1_prefix + '_acr_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    with open(opt_dir + '/temp_' + spe2_prefix + '_region.txt', 'w+') as opt:
        for eachline in store_blastedRegion_dic_in_rice_dic:
            opt.write(eachline + '\n')

    store_final_line_list = []
    with open (opt_dir + '/temp_' + spe2_prefix + '_region.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if len(col) == 3:
                store_final_line_list.append(eachline)

    with open (opt_dir + '/temp_' + spe2_prefix + '_region_modi.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_' + spe2_prefix + '_region_modi.txt > ' + opt_dir + '/temp_' + spe2_prefix + '_region_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##intersect with the Atlas
    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_' + spe2_prefix + '_region_sorted.txt' + \
          ' -b ' + ipt_rice_atlas_fl + ' > ' + opt_dir + '/temp_' + spe2_prefix + '_region_intersect_atlas_acr.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    store_acrregion_atlas_acr_celltype_dic = {}
    with open(opt_dir + '/temp_' + spe2_prefix + '_region_intersect_atlas_acr.txt', 'r') as ipt:
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


    ##collect the if the acr is root or panicle spcies
    store_acr_organ_cate_dic = {}
    all_fl_list = glob.glob(ipt_rice_all_organ_H3K27me3_dir + '/*')
    for eachfl in all_fl_list:
        mt = re.match('.+/(.+)',eachfl)
        flnm = mt.group(1)
        mt = re.match('opt_(.+)_(.+)_H3K27me3\.txt',flnm)
        organ = mt.group(2)

        with open (eachfl,'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()
                acrloc = col[0] + '_' + col[1] + '_' + col[2]
                overlapH3K27me3_cate = col[-2]
                if acrloc in store_acr_organ_cate_dic:
                    store_acr_organ_cate_dic[acrloc][organ] = overlapH3K27me3_cate
                else:
                    store_acr_organ_cate_dic[acrloc] = {}
                    store_acr_organ_cate_dic[acrloc][organ] = overlapH3K27me3_cate


    store_final_line_list = []
    with open(ipt_final_summary_fl, 'r') as ipt:
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
                        H3K27me3_organ_root_cate_list = []
                        H3K27me3_organ_panicle_cate_list = []

                        rice_atlas_acr_list = rice_atlas_acr_str.split(';')
                        for eachacr_atlas in rice_atlas_acr_list:
                            if eachacr_atlas in store_eachpeak_celltype_str_num_dic:
                                celltypenum = int(store_eachpeak_celltype_str_num_dic[eachacr_atlas]['celltypenum'])
                                celltypenum_list.append(celltypenum)


                            ##check the organ H3K27me3 cate for the atlas ACR
                            if eachacr_atlas in store_acr_organ_cate_dic:

                                if 'root' in store_acr_organ_cate_dic[eachacr_atlas]:
                                    H3K27me3_organ_root_cate = store_acr_organ_cate_dic[eachacr_atlas]['root']
                                else:
                                    H3K27me3_organ_root_cate = 'none'

                                H3K27me3_organ_root_cate_list.append(H3K27me3_organ_root_cate)

                                if 'panicle' in store_acr_organ_cate_dic[eachacr_atlas]:
                                    H3K27me3_organ_panicle_cate = store_acr_organ_cate_dic[eachacr_atlas]['panicle']
                                else:
                                    H3K27me3_organ_panicle_cate = 'none'

                                H3K27me3_organ_panicle_cate_list.append(H3K27me3_organ_panicle_cate)

                        H3K27me3_organ_root_cate_str = ','.join(H3K27me3_organ_root_cate_list)
                        H3K27me3_organ_panicle_cate_list = ','.join(H3K27me3_organ_panicle_cate_list)


                        final_celltypenum = str(int(mean(celltypenum_list)))

                        ##check the number of organs
                        store_organ_dic = {}
                        for eachacr_atlas in rice_atlas_acr_list:
                            if eachacr_atlas in store_eachpeak_celltype_str_num_dic:
                                celltype_string = store_eachpeak_celltype_str_num_dic[eachacr_atlas]['celltypestr']
                                celltype_list = celltype_string.split(',')
                                for eachcelltype in celltype_list:
                                    mt = re.match('(.+)\.+', eachcelltype)
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

                        H3K27me3_organ_root_cate_str = 'none'
                        H3K27me3_organ_panicle_cate_list = 'none'


                    final_line = eachline + '\t' + rice_atlas_acr_str + '\t' + rice_atlas_celltype_str + '\t' + organnum + '\t' + final_celltypenum + '\t' + H3K27me3_organ_root_cate_str + '\t' + H3K27me3_organ_panicle_cate_list
                    store_final_line_list.append(final_line)

                else:
                    final_line = eachline + '\t' + 'Atlas_ACR' + '\t' + 'Atlas_ACR_celltype' + '\t' + 'capture_atlas_organ_num' + '\t' + 'capture_altas_celltypenum' + '\t' + 'H3K27me3_organ_root_cate' + '\t' + 'H3K27me3_organ_panicle_cate'
                    store_final_line_list.append(final_line)

    with open(opt_dir + '/opt_' + spe1_prefix + '_to_' + spe2_prefix + '_summary_H3K27me3_add_' + spe2_prefix + '_atlas_acr.txt',
              'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



##updating 012924
##check if there are some target motif in the ACR
def subfunction_add_target_motif_in_final_summary (ipt_final_summary_fl,ipt_targret_motif_ID,ipt_spe1_motif_fl,ipt_spe2_motif_fl,opt_dir,
                                                   spe1_prefix,spe2_prefix):

    store_target_motif_loc_dic = {}
    with open (ipt_spe1_motif_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            mt = re.match('(.+)_.+',col[3])
            motifnm = mt.group(1)

            if motifnm == ipt_targret_motif_ID:

                motif_loc = col[0] + '_' + col[1] + '_' + col[2]
                store_target_motif_loc_dic[motif_loc] = 1

    with open (opt_dir + '/temp_target_motif_loc.txt','w+') as opt:
        for eachline in store_target_motif_loc_dic:
            opt.write('\t'.join(eachline.split('_')) + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_target_motif_loc.txt > ' + opt_dir + '/temp_target_motif_loc_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)


    ##for the temp target motifs, we will build a gff file loaded into the IGV
    store_final_line_list = []
    count = 0
    with open (opt_dir + '/temp_target_motif_loc_sorted.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            count += 1

            chrnm = col[0]
            st = col[1]
            ed = col[2]

            final_line = chrnm + '\t' + ipt_targret_motif_ID + '\t' + 'motif' + '\t' + st + '\t' + ed + '\t' + '.' + '\t' + '+' + '\t' + '.' + '\t' + \
                      'ID=' + ipt_targret_motif_ID + ';Name=' + ipt_targret_motif_ID + '_' + str(count)
            store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_target_motif_' + ipt_targret_motif_ID + '.gff','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    store_target_motif_loc_spe2_dic = {}
    with open (ipt_spe2_motif_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            mt = re.match('(.+)_.+', col[3])
            motifnm = mt.group(1)

            if motifnm == ipt_targret_motif_ID:

                motif_loc = col[0] + '_' + col[1] + '_' + col[2]
                store_target_motif_loc_spe2_dic[motif_loc] = 1

    with open(opt_dir + '/temp_target_motif_loc_' + spe2_prefix + '.txt', 'w+') as opt:
        for eachline in store_target_motif_loc_spe2_dic:
            opt.write('\t'.join(eachline.split('_')) + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_target_motif_loc_' + spe2_prefix + '.txt' + ' > ' + opt_dir + '/temp_target_motif_loc_' + spe2_prefix + '_sorted.txt'
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##for the temp target motifs, we will build a gff file loaded into the IGV
    store_final_line_list = []
    count = 0
    with open(opt_dir + '/temp_target_motif_loc_' + spe2_prefix + '_sorted.txt' , 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            count += 1

            chrnm = col[0]
            st = col[1]
            ed = col[2]

            final_line = chrnm + '\t' + ipt_targret_motif_ID + '\t' + 'motif' + '\t' + st + '\t' + ed + '\t' + '.' + '\t' + '+' + '\t' + '.' + '\t' + \
                         'ID=' + ipt_targret_motif_ID + ';Name=' + ipt_targret_motif_ID + '_' + str(count)
            store_final_line_list.append(final_line)

    with open(opt_dir + '/opt_target_motif_' + ipt_targret_motif_ID + '_' + spe2_prefix + '.gff', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



    ##extract the ACR information
    store_acr_loc_dic = {}
    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                acrloc = col[1]
                store_acr_loc_dic[acrloc] = 1

    with open (opt_dir + '/temp_acr_loc.txt','w+') as opt:
        for eachacrloc in store_acr_loc_dic:
            opt.write('\t'.join(eachacrloc.split('_')) + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + opt_dir + '/temp_acr_loc.txt > ' + opt_dir + '/temp_acr_loc_sorted.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    cmd = 'bedtools intersect -wa -wb -a ' + opt_dir + '/temp_target_motif_loc_sorted.txt' + \
          ' -b ' + opt_dir + '/temp_acr_loc_sorted.txt > ' + opt_dir + '/temp_intersect_motif_acr.txt'
    print(cmd)
    subprocess.call(cmd,shell=True)

    store_acr_loc_include_target_motif_dic = {}
    with open (opt_dir + '/temp_intersect_motif_acr.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            motifloc = col[0] + '_' + col[1] + '_' + col[2]
            acrloc = col[3] + '_' + col[4] + '_' + col[5]
            store_acr_loc_include_target_motif_dic[acrloc] = 1

    store_final_line_list = []
    count = 0
    with open (ipt_final_summary_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            acrloc = col[1]

            count += 1
            if count != 1:
                if acrloc in store_acr_loc_include_target_motif_dic:
                    motif_cover = 'yes'
                else:
                    motif_cover = 'no'

                final_line = eachline + '\t' + motif_cover
                store_final_line_list.append(final_line)

            else:
                final_line = eachline + '\t' + 'motif_cover'
                store_final_line_list.append(final_line)

    with open (opt_dir + '/opt_final_summary_file_add_' + spe1_prefix + '_' + spe2_prefix + '.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



##updating 013124
##we will check the 





















