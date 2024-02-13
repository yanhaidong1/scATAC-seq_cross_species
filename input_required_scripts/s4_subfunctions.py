#!/usr/bin/env python

##updating
##updating 112123 we will find the shared acr not consuder the H3K27me3

import subprocess





def subfunction_iden_shared_acr_under_H3K27me3 (ipt_spe1_to_spe2_final_fl,ipt_spe1_to_spe3_final_fl,
                                                ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,ipt_spe3_acr_cate_fl,opt_dir,
                                                ipt_spe1_prefix,ipt_spe2_prefix,ipt_spe3_prefix):

    store_spe1_acr_cate_dic = {}
    with open (ipt_spe1_acr_cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_spe1_acr_cate_dic[col[0]] = col[1]

    store_spe2_acr_cate_dic = {}
    with open (ipt_spe2_acr_cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_spe2_acr_cate_dic[col[0]] = col[1]

    store_spe3_acr_cate_dic = {}
    with open (ipt_spe3_acr_cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_spe3_acr_cate_dic[col[0]] = col[1]

    store_spe1_spe2_acr_cate_line = []

    store_spe1_acr_to_spe2_acr_BroadH3K27me3_dic = {}
    store_spe1_acr_to_spe2_acr_BroadToRestrictedH3K27me3_dic = {}
    count = 0
    with open (ipt_spe1_to_spe2_final_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split()
                spe1_acr = col[1]
                spe2_acr = col[5]
                spe1_H3K27me3cate = col[2]
                spe2_H3K27me3cate = col[6]

                syntenic_regionID = col[11]
                spe1_gene_direction = col[10]
                spe2_gene_direction = col[13]

                if len(syntenic_regionID.split(',')) == 1:

                    if spe1_gene_direction == 'OneDirection' and spe2_gene_direction == 'OneDirection':

                        if spe1_H3K27me3cate == 'BroadInflankH3K27me3peak' and spe2_H3K27me3cate == 'BroadInflankH3K27me3peak':

                            if spe1_acr in store_spe1_acr_cate_dic and spe2_acr in store_spe2_acr_cate_dic:
                                if store_spe1_acr_cate_dic[spe1_acr] == 'Distal' or store_spe1_acr_cate_dic[spe1_acr] == 'Proximal':
                                    if store_spe2_acr_cate_dic[spe2_acr] == 'Distal' or store_spe2_acr_cate_dic[spe2_acr] == 'Proximal':

                                        if spe1_acr in store_spe1_acr_to_spe2_acr_BroadH3K27me3_dic:
                                            store_spe1_acr_to_spe2_acr_BroadH3K27me3_dic[spe1_acr][spe2_acr] = syntenic_regionID
                                        else:
                                            store_spe1_acr_to_spe2_acr_BroadH3K27me3_dic[spe1_acr] = {}
                                            store_spe1_acr_to_spe2_acr_BroadH3K27me3_dic[spe1_acr][spe2_acr] = syntenic_regionID


                        if spe1_H3K27me3cate == 'BroadInflankH3K27me3peak' and spe2_H3K27me3cate == 'RestrictedInflankH3K27me3peak':

                            if spe1_acr in store_spe1_acr_cate_dic and spe2_acr in store_spe2_acr_cate_dic:
                                if store_spe1_acr_cate_dic[spe1_acr] == 'Distal' or store_spe1_acr_cate_dic[spe1_acr] == 'Proximal':
                                    if store_spe2_acr_cate_dic[spe2_acr] == 'Distal' or store_spe2_acr_cate_dic[spe2_acr] == 'Proximal':

                                        if spe1_acr in store_spe1_acr_to_spe2_acr_BroadToRestrictedH3K27me3_dic:
                                            store_spe1_acr_to_spe2_acr_BroadToRestrictedH3K27me3_dic[spe1_acr][spe2_acr] = syntenic_regionID
                                        else:
                                            store_spe1_acr_to_spe2_acr_BroadToRestrictedH3K27me3_dic[spe1_acr] = {}
                                            store_spe1_acr_to_spe2_acr_BroadToRestrictedH3K27me3_dic[spe1_acr][spe2_acr] = syntenic_regionID



                if spe1_acr in store_spe1_acr_cate_dic and spe2_acr in store_spe2_acr_cate_dic:
                    spe1_cate = store_spe1_acr_cate_dic[spe1_acr]
                    spe2_cate = store_spe2_acr_cate_dic[spe2_acr]

                    final_line = spe1_acr + '\t' + spe1_cate + '\t' + spe1_H3K27me3cate + '\t' + spe2_acr + '\t' + spe2_cate + '\t' + spe2_H3K27me3cate
                    store_spe1_spe2_acr_cate_line.append(final_line)

    print(store_spe1_acr_to_spe2_acr_BroadToRestrictedH3K27me3_dic)

    with open (opt_dir + '/opt_' + ipt_spe1_prefix + '_' + ipt_spe2_prefix + '_acr_gene_cate.txt','w+') as opt:
        for eachline in store_spe1_spe2_acr_cate_line:
            opt.write(eachline + '\n')

    store_spe1_spe3_acr_cate_line = []
    store_spe1_acr_to_spe3_acr_BroadH3K27me3_dic = {}
    store_spe1_acr_to_spe3_acr_BroadToRestrictedH3K27me3_dic = {}
    count = 0
    with open (ipt_spe1_to_spe3_final_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split()
                spe1_acr = col[1]
                spe3_acr = col[5]
                spe1_H3K27me3cate = col[2]
                spe3_H3K27me3cate = col[6]
                syntenic_regionID = col[11]
                spe1_gene_direction = col[10]
                spe2_gene_direction = col[13]

                if len(syntenic_regionID.split(',')) == 1:

                    if spe1_gene_direction == 'OneDirection' and spe2_gene_direction == 'OneDirection':

                        if spe1_H3K27me3cate == 'BroadInflankH3K27me3peak' and spe3_H3K27me3cate == 'BroadInflankH3K27me3peak':

                            if spe1_acr in store_spe1_acr_cate_dic and spe3_acr in store_spe3_acr_cate_dic:
                                if store_spe1_acr_cate_dic[spe1_acr] == 'Distal' or store_spe1_acr_cate_dic[spe1_acr] == 'Proximal':
                                    if store_spe3_acr_cate_dic[spe3_acr] == 'Distal' or store_spe3_acr_cate_dic[spe3_acr] == 'Proximal':

                                        if spe1_acr in store_spe1_acr_to_spe3_acr_BroadH3K27me3_dic:
                                            store_spe1_acr_to_spe3_acr_BroadH3K27me3_dic[spe1_acr][spe3_acr] = syntenic_regionID
                                        else:
                                            store_spe1_acr_to_spe3_acr_BroadH3K27me3_dic[spe1_acr] = {}
                                            store_spe1_acr_to_spe3_acr_BroadH3K27me3_dic[spe1_acr][spe3_acr] = syntenic_regionID


                        if spe1_H3K27me3cate == 'BroadInflankH3K27me3peak' and spe3_H3K27me3cate == 'RestrictedInflankH3K27me3peak':

                            if spe1_acr in store_spe1_acr_cate_dic and spe3_acr in store_spe3_acr_cate_dic:
                                if store_spe1_acr_cate_dic[spe1_acr] == 'Distal' or store_spe1_acr_cate_dic[
                                    spe1_acr] == 'Proximal':
                                    if store_spe3_acr_cate_dic[spe3_acr] == 'Distal' or store_spe3_acr_cate_dic[
                                        spe3_acr] == 'Proximal':

                                        if spe1_acr in store_spe1_acr_to_spe3_acr_BroadToRestrictedH3K27me3_dic:
                                            store_spe1_acr_to_spe3_acr_BroadToRestrictedH3K27me3_dic[spe1_acr][
                                                spe3_acr] = syntenic_regionID
                                        else:
                                            store_spe1_acr_to_spe3_acr_BroadToRestrictedH3K27me3_dic[spe1_acr] = {}
                                            store_spe1_acr_to_spe3_acr_BroadToRestrictedH3K27me3_dic[spe1_acr][
                                                spe3_acr] = syntenic_regionID


                if spe1_acr in store_spe1_acr_cate_dic and spe3_acr in store_spe3_acr_cate_dic:
                    spe1_cate = store_spe1_acr_cate_dic[spe1_acr]
                    spe3_cate = store_spe3_acr_cate_dic[spe3_acr]

                    final_line = spe1_acr + '\t' + spe1_cate + '\t' + spe1_H3K27me3cate + '\t' + spe3_acr + '\t' + spe3_cate + '\t' + spe3_H3K27me3cate
                    store_spe1_spe3_acr_cate_line.append(final_line)

    print(store_spe1_acr_to_spe3_acr_BroadToRestrictedH3K27me3_dic)

    with open (opt_dir + '/opt_' + ipt_spe1_prefix + '_' + ipt_spe3_prefix + '_acr_gene_cate.txt','w+') as opt:
        for eachline in store_spe1_spe3_acr_cate_line:
            opt.write(eachline + '\n')


    ##check how many of them could be identified in both species
    store_final_overlapped_acr_line_BroadToBroad_H3K27me3_list = []
    for eachspe1_acr in store_spe1_acr_to_spe2_acr_BroadH3K27me3_dic:
        if eachspe1_acr in store_spe1_acr_to_spe3_acr_BroadH3K27me3_dic:

            spe2_acr_dic = store_spe1_acr_to_spe2_acr_BroadH3K27me3_dic[eachspe1_acr]
            spe3_acr_dic = store_spe1_acr_to_spe3_acr_BroadH3K27me3_dic[eachspe1_acr]

            spe2_acr_str = ','.join(list(spe2_acr_dic.keys()))
            spe3_acr_str = ','.join(list(spe3_acr_dic.keys()))

            spe2_region_list = []
            for eachacr in spe2_acr_dic:
                spe2_region_list.append(spe2_acr_dic[eachacr])
            spe2_region_str = ','.join(spe2_region_list)

            spe3_region_list = []
            for eachacr in spe3_acr_dic:
                spe3_region_list.append(spe3_acr_dic[eachacr])
            spe3_region_str = ','.join(spe3_region_list)

            final_line = eachspe1_acr + '\t' + spe2_acr_str + '\t' + spe2_region_str + '\t' + spe3_acr_str + '\t' + spe3_region_str
            store_final_overlapped_acr_line_BroadToBroad_H3K27me3_list.append(final_line)


    store_final_overlapped_acr_line_BroadToRestricted_H3K27me3_list = []
    for eachspe1_acr in store_spe1_acr_to_spe2_acr_BroadToRestrictedH3K27me3_dic:
        if eachspe1_acr in store_spe1_acr_to_spe3_acr_BroadToRestrictedH3K27me3_dic:

            spe2_acr_dic = store_spe1_acr_to_spe2_acr_BroadToRestrictedH3K27me3_dic[eachspe1_acr]
            spe3_acr_dic = store_spe1_acr_to_spe3_acr_BroadToRestrictedH3K27me3_dic[eachspe1_acr]

            spe2_acr_str = ','.join(list(spe2_acr_dic.keys()))
            spe3_acr_str = ','.join(list(spe3_acr_dic.keys()))

            spe2_region_list = []
            for eachacr in spe2_acr_dic:
                spe2_region_list.append(spe2_acr_dic[eachacr])
            spe2_region_str = ','.join(spe2_region_list)

            spe3_region_list = []
            for eachacr in spe3_acr_dic:
                spe3_region_list.append(spe3_acr_dic[eachacr])
            spe3_region_str = ','.join(spe3_region_list)

            final_line = eachspe1_acr + '\t' + spe2_acr_str + '\t' + spe2_region_str + '\t' + spe3_acr_str + '\t' + spe3_region_str
            store_final_overlapped_acr_line_BroadToRestricted_H3K27me3_list.append(final_line)



    return (store_final_overlapped_acr_line_BroadToBroad_H3K27me3_list,store_final_overlapped_acr_line_BroadToRestricted_H3K27me3_list)




################################################
def store_spe_acr_cate_dic (ipt_spe_acr_cate_fl):

    store_spe1_acr_cate_dic = {}
    with open (ipt_spe_acr_cate_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_spe1_acr_cate_dic[col[0]] = col[1]

    return (store_spe1_acr_cate_dic)

def store_spe1_acr_to_spe2 (ipt_spe1_to_spe2_final_fl,
                            store_spe1_acr_cate_dic,
                            store_spe2_acr_cate_dic,
                            ipt_spe1_prefix,ipt_spe2_prefix,opt_dir):

    store_spe1_spe2_acr_cate_line = []
    store_spe1_acr_to_spe2_acr_dic = {}
    count = 0
    with open(ipt_spe1_to_spe2_final_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split()
                spe1_acr = col[1]
                spe2_acr = col[5]

                syntenic_regionID = col[11]
                spe1_gene_direction = col[10]
                spe2_gene_direction = col[13]

                spe1_celltype = col[8]
                spe2_celltype = col[14]


                if len(syntenic_regionID.split(',')) == 1:

                    if spe1_gene_direction == 'OneDirection' and spe2_gene_direction == 'OneDirection':

                        if spe1_celltype == 'broadly_accessible' and spe2_celltype == 'broadly_accessible':

                            if spe1_acr in store_spe1_acr_cate_dic and spe2_acr in store_spe2_acr_cate_dic:
                                if store_spe1_acr_cate_dic[spe1_acr] == 'Distal' or store_spe1_acr_cate_dic[
                                    spe1_acr] == 'Proximal':
                                    if store_spe2_acr_cate_dic[spe2_acr] == 'Distal' or store_spe2_acr_cate_dic[
                                        spe2_acr] == 'Proximal':

                                        if spe1_acr in store_spe1_acr_to_spe2_acr_dic:
                                            store_spe1_acr_to_spe2_acr_dic[spe1_acr][
                                                spe2_acr] = syntenic_regionID
                                        else:
                                            store_spe1_acr_to_spe2_acr_dic[spe1_acr] = {}
                                            store_spe1_acr_to_spe2_acr_dic[spe1_acr][
                                                spe2_acr] = syntenic_regionID

                if spe1_acr in store_spe1_acr_cate_dic and spe2_acr in store_spe2_acr_cate_dic:
                    spe1_cate = store_spe1_acr_cate_dic[spe1_acr]
                    spe2_cate = store_spe2_acr_cate_dic[spe2_acr]

                    final_line = spe1_acr + '\t' + spe1_cate  + '\t' + spe2_acr + '\t' + spe2_cate
                    store_spe1_spe2_acr_cate_line.append(final_line)

    with open (opt_dir + '/opt_' + ipt_spe1_prefix + '_' + ipt_spe2_prefix + '_acr_gene_cate.txt','w+') as opt:
        for eachline in store_spe1_spe2_acr_cate_line:
            opt.write(eachline + '\n')

    return (store_spe1_acr_to_spe2_acr_dic)


def subfunction_iden_shared_acr_five_species (ipt_spe1_to_spe2_final_fl,ipt_spe1_to_spe3_final_fl,ipt_spe1_to_spe4_final_fl,ipt_spe1_to_spe5_final_fl,
                                                ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,ipt_spe3_acr_cate_fl,ipt_spe4_acr_cate_fl,ipt_spe5_acr_cate_fl,
                                                opt_dir,
                                                ipt_spe1_prefix,ipt_spe2_prefix,ipt_spe3_prefix,ipt_spe4_prefix,ipt_spe5_prefix):

    ##store the spe acr cate dic
    store_spe1_acr_cate_dic = store_spe_acr_cate_dic(ipt_spe1_acr_cate_fl)
    store_spe2_acr_cate_dic = store_spe_acr_cate_dic(ipt_spe2_acr_cate_fl)
    store_spe3_acr_cate_dic = store_spe_acr_cate_dic(ipt_spe3_acr_cate_fl)
    store_spe4_acr_cate_dic = store_spe_acr_cate_dic(ipt_spe4_acr_cate_fl)
    store_spe5_acr_cate_dic = store_spe_acr_cate_dic(ipt_spe5_acr_cate_fl)

    store_spe1_acr_to_spe2_acr_dic = store_spe1_acr_to_spe2(ipt_spe1_to_spe2_final_fl,
                                                                   store_spe1_acr_cate_dic,
                                                                   store_spe2_acr_cate_dic,
                                                                   ipt_spe1_prefix,ipt_spe2_prefix,opt_dir)

    store_spe1_acr_to_spe3_acr_dic = store_spe1_acr_to_spe2(ipt_spe1_to_spe3_final_fl,
                                                                   store_spe1_acr_cate_dic,
                                                                   store_spe3_acr_cate_dic,
                                                            ipt_spe1_prefix, ipt_spe3_prefix,opt_dir)

    store_spe1_acr_to_spe4_acr_dic = store_spe1_acr_to_spe2(ipt_spe1_to_spe4_final_fl,
                                                                   store_spe1_acr_cate_dic,
                                                                   store_spe4_acr_cate_dic,
                                                            ipt_spe1_prefix, ipt_spe4_prefix,opt_dir)

    store_spe1_acr_to_spe5_acr_dic = store_spe1_acr_to_spe2(ipt_spe1_to_spe5_final_fl,
                                                                   store_spe1_acr_cate_dic,
                                                                   store_spe5_acr_cate_dic,
                                                            ipt_spe1_prefix, ipt_spe5_prefix,opt_dir)

    ##check how many of them could be identified in both species
    store_final_overlapped_acr_line_five_species_list = []
    for eachspe1_acr in store_spe1_acr_to_spe2_acr_dic:
        if eachspe1_acr in store_spe1_acr_to_spe3_acr_dic:
            if eachspe1_acr in store_spe1_acr_to_spe4_acr_dic:
                if eachspe1_acr in store_spe1_acr_to_spe5_acr_dic:

                    spe2_acr_dic = store_spe1_acr_to_spe2_acr_dic[eachspe1_acr]
                    spe3_acr_dic = store_spe1_acr_to_spe3_acr_dic[eachspe1_acr]
                    spe4_acr_dic = store_spe1_acr_to_spe4_acr_dic[eachspe1_acr]
                    spe5_acr_dic = store_spe1_acr_to_spe5_acr_dic[eachspe1_acr]

                    spe2_acr_str = ','.join(list(spe2_acr_dic.keys()))
                    spe3_acr_str = ','.join(list(spe3_acr_dic.keys()))
                    spe4_acr_str = ','.join(list(spe4_acr_dic.keys()))
                    spe5_acr_str = ','.join(list(spe5_acr_dic.keys()))

                    spe2_region_list = []
                    for eachacr in spe2_acr_dic:
                        spe2_region_list.append(spe2_acr_dic[eachacr])
                    spe2_region_str = ','.join(spe2_region_list)

                    spe3_region_list = []
                    for eachacr in spe3_acr_dic:
                        spe3_region_list.append(spe3_acr_dic[eachacr])
                    spe3_region_str = ','.join(spe3_region_list)

                    spe4_region_list = []
                    for eachacr in spe4_acr_dic:
                        spe4_region_list.append(spe4_acr_dic[eachacr])
                    spe4_region_str = ','.join(spe4_region_list)

                    spe5_region_list = []
                    for eachacr in spe5_acr_dic:
                        spe5_region_list.append(spe5_acr_dic[eachacr])
                    spe5_region_str = ','.join(spe5_region_list)


                    final_line = eachspe1_acr + '\t' + spe2_acr_str + '\t' + spe2_region_str + '\t' + spe3_acr_str + '\t' + spe3_region_str + \
                                '\t' + spe4_acr_str + '\t' + spe4_region_str + '\t' + spe5_acr_str + '\t' + spe5_region_str

                    store_final_overlapped_acr_line_five_species_list.append(final_line)


    return (store_final_overlapped_acr_line_five_species_list)


################################################
def store_spe1_acr_to_spe2_celltype_specific (ipt_spe1_to_spe2_final_fl,
                            store_spe1_acr_cate_dic,
                            store_spe2_acr_cate_dic,
                            ipt_spe1_prefix,ipt_spe2_prefix,ipt_target_celltype_nm,opt_dir):

    store_spe1_spe2_acr_cate_line = []
    store_spe1_acr_to_spe2_acr_dic = {}
    count = 0
    with open(ipt_spe1_to_spe2_final_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split()
                spe1_acr = col[1]
                spe2_acr = col[5]

                syntenic_regionID = col[11]
                spe1_gene_direction = col[10]
                spe2_gene_direction = col[13]

                spe1_celltype = col[8]
                spe2_celltype = col[14]


                if len(syntenic_regionID.split(',')) == 1:

                    if spe1_gene_direction == 'OneDirection' and spe2_gene_direction == 'OneDirection':

                        if ipt_target_celltype_nm in spe1_celltype and ipt_target_celltype_nm in spe2_celltype:

                            if spe1_acr in store_spe1_acr_cate_dic and spe2_acr in store_spe2_acr_cate_dic:
                                if store_spe1_acr_cate_dic[spe1_acr] == 'Distal' or store_spe1_acr_cate_dic[
                                    spe1_acr] == 'Proximal':
                                    if store_spe2_acr_cate_dic[spe2_acr] == 'Distal' or store_spe2_acr_cate_dic[
                                        spe2_acr] == 'Proximal':

                                        if spe1_acr in store_spe1_acr_to_spe2_acr_dic:
                                            store_spe1_acr_to_spe2_acr_dic[spe1_acr][
                                                spe2_acr] = syntenic_regionID
                                        else:
                                            store_spe1_acr_to_spe2_acr_dic[spe1_acr] = {}
                                            store_spe1_acr_to_spe2_acr_dic[spe1_acr][
                                                spe2_acr] = syntenic_regionID

                if spe1_acr in store_spe1_acr_cate_dic and spe2_acr in store_spe2_acr_cate_dic:
                    spe1_cate = store_spe1_acr_cate_dic[spe1_acr]
                    spe2_cate = store_spe2_acr_cate_dic[spe2_acr]

                    final_line = spe1_acr + '\t' + spe1_cate  + '\t' + spe2_acr + '\t' + spe2_cate
                    store_spe1_spe2_acr_cate_line.append(final_line)

    with open (opt_dir + '/opt_' + ipt_spe1_prefix + '_' + ipt_spe2_prefix + '_acr_gene_cate.txt','w+') as opt:
        for eachline in store_spe1_spe2_acr_cate_line:
            opt.write(eachline + '\n')

    return (store_spe1_acr_to_spe2_acr_dic)


def subfunction_iden_celltype_specific_ACR_correspond (ipt_spe1_to_spe2_final_fl,
                                                       ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,
                                                       ipt_spe1_prefix,ipt_spe2_prefix,
                                                       ipt_target_celltype_nm,opt_dir):

    store_spe1_acr_cate_dic = store_spe_acr_cate_dic(ipt_spe1_acr_cate_fl)
    store_spe2_acr_cate_dic = store_spe_acr_cate_dic(ipt_spe2_acr_cate_fl)

    store_spe1_acr_to_spe2_acr_dic = store_spe1_acr_to_spe2_celltype_specific(ipt_spe1_to_spe2_final_fl,
                                             store_spe1_acr_cate_dic,
                                             store_spe2_acr_cate_dic,
                                             ipt_spe1_prefix, ipt_spe2_prefix, ipt_target_celltype_nm, opt_dir)

    ##check how many of them could be identified in both species
    store_final_overlapped_acr_line_celltype_list = []
    for eachspe1_acr in store_spe1_acr_to_spe2_acr_dic:

        spe2_acr_dic = store_spe1_acr_to_spe2_acr_dic[eachspe1_acr]

        spe2_acr_str = ','.join(list(spe2_acr_dic.keys()))


        spe2_region_list = []
        for eachacr in spe2_acr_dic:
            spe2_region_list.append(spe2_acr_dic[eachacr])
        spe2_region_str = ','.join(spe2_region_list)

        final_line = eachspe1_acr + '\t' + spe2_acr_str + '\t' + spe2_region_str + '\t' + ipt_target_celltype_nm

        store_final_overlapped_acr_line_celltype_list.append(final_line)

    return (store_final_overlapped_acr_line_celltype_list)




















