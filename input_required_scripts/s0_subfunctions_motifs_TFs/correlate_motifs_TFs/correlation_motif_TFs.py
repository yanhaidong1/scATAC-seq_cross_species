#!/usr/bin/env python

##updating 112723 we will change to a version that has already know the corresponding TFs
##updating 053123 we will allow the positive correlation to be higher than 0.1
##updating 031323 we will add a new function to generate another correlation and consider each ortholog
##updating 122322 set another case to show the corresponding
##updating 122222
##1) generate a new PerM file with correspnding TF name for all the motifs
##2) filter TF to generate a new matrix that corresponds to the motifs
##3) we will use the R to generate the correlation


##this script is to correlate motif enrichment FC and TFs
##this correlation is based on the cluster other than the cell since the fimoModel reg enrichment analysis cannot obtain the acc of each cell

##since one At TF may correspond many different rice genes, we need to check the correlation value for each of them and select the largest one

import re
import glob
import sys
import subprocess
import os
#from scipy.stats.stats import spearmanr
from scipy.stats.stats import pearsonr

input_blasttf_fl = sys.argv[1] ##opt_flt03_addtffam_fltdtf_blast.txt
##

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/11_TFmotif_analysis_040522/output_dir_072522/step02_motif_TF_corresponding_fimoModel/s1_make_blast_TFs_dir/opt_flt015_addtffam_fltdtf_blast.txt


input_Atgene_TFid_fl = sys.argv[2]



input_Atgene_motif_fl = sys.argv[3]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/11_TFmotif_analysis_040522/output_dir_072522/step02_motif_TF_corresponding_fimoModel/s2_correspond_motif_gene_dir/opt_Atgene_motif.txt

input_norm_gene_act_fl = sys.argv[4]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/01_summary_of_QC_020322/09_3_get_PerM_allorgans_from_raw_gene_count_sparse_120622/output_dir_downsampling_70kcells_122122/opt_allorgans_perM_genes_accessible_clusters.txt

##/scratch/hy17471/rice_altlas_scATAC_seq_042021/01_summary_of_QC_020322/09_2_get_PerM_allorgans_normalized_acc_081322/output_dir/opt_allorgans_perM_genes_accessible_clusters.txt

input_norm_motif_act_fl = sys.argv[5]
##this is the average motif acr

#input_motif_deviation_fl = sys.argv[4]

input_corr_R_script_fl = sys.argv[6]


input_open_compare_corr_rate = sys.argv[7]
##if we open this function, we will set if the corr in coefficient is lower than 0.2, we will select the largest from negative correlation.


##input_motif_fc_enrichment_score_fl = sys.argv[4]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/11_TFmotif_analysis_040522/output_dir_072522/step01_celltype_enrich_motif_dir/s7_do_enrichment_dir_opt2/opt_motif_enrichment_cluster.txt


##/scratch/hy17471/soybean_scATAC_100120/pipeline_analysis_110220/add_11_TF_motif_cor_011221/02_add_jasprmotif_to_blast_012521/output_dir_012921/opt_Atgene_motif.txt
#input_norm_gene_act_fl = sys.argv[3]
## updating if we obtain new normalized act sparse
##/scratch/hy17471/soybean_scATAC_100120/pipeline_analysis_110220/add_06_geneAccessibility_UMAP_analysis_regmodel2_res01_shift_011421/all.normalizedActivity.sparse
##we do not consider the binary since the previous UMAP utilizes the non-binary to have an analysis
#input_motif_dev_score_fl = sys.argv[4]
##directly generate sparse
##/scratch/hy17471/soybean_scATAC_100120/pipeline_analysis_110220/add_11_chromVAR_motif_dev_010221/02_run_chromVAR_012921/opt_motif_deviations.sparse

input_output_dir = sys.argv[8]



def corr_TFacc_motifdev (input_blasttf_fl,input_Atgene_TFid_fl,input_Atgene_motif_fl,input_norm_gene_act_fl,input_norm_motif_act_fl,input_output_dir,
                         input_corr_R_script_fl,input_open_compare_corr_rate):

    ##store the Atgene TFid dic
    store_TFid_Atgene_dic = {}
    with open (input_Atgene_TFid_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            store_TFid_Atgene_dic[col[4]] = col[3]


    ##find the TFs needs to be analyzed
    store_atgene_withmotif_infor_dic = {}
    with open (input_Atgene_motif_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            gene_nm = col[0]
            motif_nm = col[1] ##also contains MA1224.1_CBF1;MA1669.1_DREB1B
            upper_gene_nm = gene_nm.upper()
            store_atgene_withmotif_infor_dic[upper_gene_nm] = motif_nm

    ##one motif may corrrespond to multipe atgenes
    store_motif_atgene_dic = {}
    with open(input_Atgene_motif_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            atgene = col[0]
            motif_nm = col[1]  ##also contains MA1224.1_CBF1;MA1669.1_DREB1B

            motif_nm_list = motif_nm.split(';')

            for eachmotif in motif_nm_list:
                if eachmotif in store_motif_atgene_dic:
                    store_motif_atgene_dic[eachmotif][atgene] = 1
                else:
                    store_motif_atgene_dic[eachmotif] = {}
                    store_motif_atgene_dic[eachmotif][atgene] = 1

    ##one atgene may correspond to multiple rice gene
    store_atgene_wmotif_soylist_dic = {} ##key is the at gene and value is a dic to store the soybean gene
    count = 0
    with open (input_blasttf_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
            #if not eachline.startswith('chr'):

                col = eachline.strip().split()
                if re.match('.+\.MSUv7\.0',col[3]):
                    mt = re.match('(.+)\.MSUv7\.0',col[3])
                    spe_gene_ID = mt.group(1)
                else:
                    spe_gene_ID = col[3]


                At_gene_list = []
                if '__' in col[4]:

                    all_ortho_col = col[4].split('__')
                    for eachortho in all_ortho_col:

                        if re.match('(.+)_ortho.+',eachortho):

                            mt = re.match('(.+)_ortho.+',eachortho)
                            TFid = mt.group(1)

                        else:
                            TFid = eachortho

                        Atgene = store_TFid_Atgene_dic[TFid]
                        At_gene_list.append(Atgene)

                else:
                    if '_ortho' in col[4]:
                        mt = re.match('(.+)_ortho.+', col[4])
                        TFid = mt.group(1)
                        Atgene = store_TFid_Atgene_dic[TFid]
                        At_gene_list.append(Atgene)

                    else:
                        TFid = col[4]
                        Atgene = store_TFid_Atgene_dic[TFid]
                        At_gene_list.append(Atgene)

                for eachatgene in At_gene_list:

                    if eachatgene in store_atgene_wmotif_soylist_dic:
                        store_atgene_wmotif_soylist_dic[eachatgene][spe_gene_ID] = 1
                    else:
                        store_atgene_wmotif_soylist_dic[eachatgene] = {}
                        store_atgene_wmotif_soylist_dic[eachatgene][spe_gene_ID] = 1


    ##write out target soybean genes
    store_target_soygene_dic = {}
    for eachgene in store_atgene_wmotif_soylist_dic:
        for eachsoygene in store_atgene_wmotif_soylist_dic[eachgene]:
            store_target_soygene_dic[eachsoygene] = 1
    with open (input_output_dir + '/temp_target_QueryGene_list.txt','w+') as opt:
        for eachgene in store_target_soygene_dic:
            opt.write(eachgene + '\n')

    ##now we will corresponding the TF and motif names
    ##we first corresponding the motif name to the TF
    store_corresponding_motif_ricegene_line_list = []
    count = 0
    with open (input_norm_motif_act_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1

            if count != 1:
                col = eachline.strip().split()
                mt = re.match('(.+)_.+',col[0])
                motifnm = mt.group(1)

                if motifnm in store_motif_atgene_dic:

                    atgene_dic = store_motif_atgene_dic[motifnm]

                    for eachatgene in atgene_dic:

                        if eachatgene in store_atgene_wmotif_soylist_dic:
                            riccegene_dic = store_atgene_wmotif_soylist_dic[eachatgene]

                            final_line = motifnm + '\t' + ';'.join(list(riccegene_dic.keys()))
                            store_corresponding_motif_ricegene_line_list.append(final_line)

    with open (input_output_dir + '/temp_corresponding_motif_ricegene.txt','w+') as opt:
        for eachline in store_corresponding_motif_ricegene_line_list:
            opt.write(eachline + '\n')

    ##udpating 060223
    ##we will filter the genes in the gene act fl based on the target total genes
    store_correspond_motif_gene_list_dic = {}
    with open (input_output_dir + '/temp_corresponding_motif_ricegene.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if len(col) == 2:
                all_gene_list = col[1].split(';')
                for eachgene in all_gene_list:
                    store_correspond_motif_gene_list_dic[eachgene] = 1

    with open (input_output_dir + '/temp_gene_correspond_to_motif.txt','w+') as opt:
        for eachgene in store_correspond_motif_gene_list_dic:
            opt.write(eachgene + '\n')


    store_final_line_list = []
    count = 0
    with open (input_norm_gene_act_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split()
                geneID = col[0]
                if geneID in store_correspond_motif_gene_list_dic:
                    store_final_line_list.append(eachline)
            else:
                store_final_line_list.append(eachline)

    with open (input_output_dir + '/opt_TF_average_smooth_accessible_clusters_fltgenes.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##do the correlation between two matrix
    ##the basic ideas is to select the homolog with the greatest Pearson correlation coefficient with respect to the motif devaiation score.
    cmd = 'Rscript ' + input_corr_R_script_fl + \
          ' ' + input_output_dir + '/opt_TF_average_smooth_accessible_clusters_fltgenes.txt' + \
          ' ' + input_norm_motif_act_fl + \
          ' ' + input_output_dir
    print(cmd)
    subprocess.call(cmd,shell=True)

    opt_corr_dt = input_output_dir + '/opt_corr_dt.txt'

    store_tfgene_motif_dic = {}
    with open (input_output_dir + '/temp_corresponding_motif_ricegene.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            geneid = col[1]
            genelist = geneid.split(';')
            for eachgene in genelist:
                if eachgene in store_tfgene_motif_dic:
                    store_tfgene_motif_dic[eachgene][col[0]] = 1
                else:
                    store_tfgene_motif_dic[eachgene] = {}
                    store_tfgene_motif_dic[eachgene][col[0]] = 1

    store_gene_motif_line_dic = {}
    store_motif_order_dic = {}
    count = 0
    with open(opt_corr_dt, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                store_gene_motif_line_dic[col[0]] = eachline
            else:
                for i in range(len(col)):
                    mt = re.match('(.+)_.+',col[i])
                    motifnm = mt.group(1)
                    store_motif_order_dic[motifnm] = i + 1

    ##find the most match motif name
    store_final_line_list = []
    for eachgene in store_tfgene_motif_dic:
        motif_dic = store_tfgene_motif_dic[eachgene]
        motif_str = ','.join(list(motif_dic.keys()))

        ##now we will extract information from opt_corr_dt
        if len(list(motif_dic.keys())) != 1:

            if eachgene in store_gene_motif_line_dic:

                motif_line = store_gene_motif_line_dic[eachgene]
                col = motif_line.strip().split()

                store_motif_val_dic = {}
                motif_val_list = []
                for eachmotif in motif_dic:

                    motif_value_order = store_motif_order_dic[eachmotif]
                    motif_val = col[motif_value_order]

                    store_motif_val_dic[eachmotif] = float(motif_val)
                    motif_val_list.append(motif_val)

                ##updating 122422
                if input_open_compare_corr_rate == 'no':

                    max_motif = max(store_motif_val_dic, key=store_motif_val_dic.get)
                    final_line = eachgene + '\t' + max_motif + '\t' + motif_str + '\t' + ','.join(motif_val_list)

                else:

                    temp_max_motif = max(store_motif_val_dic, key=store_motif_val_dic.get)

                    if store_motif_val_dic[temp_max_motif] >= 0.1:
                        final_line = eachgene + '\t' + temp_max_motif + '\t' + motif_str + '\t' + ','.join(motif_val_list)
                    else:
                        ##we will select the max of negative ones
                        store_abs_motif_val_dic = {}
                        for eachabs_motif in store_motif_val_dic:
                            store_abs_motif_val_dic[eachabs_motif] = abs(store_motif_val_dic[eachabs_motif])
                        temp_abs_max_motif = max(store_abs_motif_val_dic, key=store_abs_motif_val_dic.get)

                        final_line = eachgene + '\t' + temp_abs_max_motif + '\t' + motif_str + '\t' + ','.join(motif_val_list)

                store_final_line_list.append(final_line)

            else:
                print('Warning: Gene ' + eachgene + ' not in store_gene_motif_line_dic')

        else:
            final_line = eachgene + '\t' + motif_str + '\t' + motif_str + '\t' + 'na'
            store_final_line_list.append(final_line)

    with open (input_output_dir + '/temp_corresponding_motif_ricegene_flt1.txt','w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')


    ##since there are cases one motif may have two TFs
    ##we need to find the best match
    opt_corr_dt = input_output_dir + '/opt_corr_dt_rev.txt'

    store_motif_gene_dic = {}
    with open(input_output_dir + '/temp_corresponding_motif_ricegene_flt1.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            geneid = col[0]
            motif = col[1]

            if motif in store_motif_gene_dic:
                store_motif_gene_dic[motif][geneid] = 1
            else:
                store_motif_gene_dic[motif] = {}
                store_motif_gene_dic[motif][geneid] = 1


    store_motif_gene_line_dic = {}
    store_gene_order_dic = {}
    count = 0
    with open(opt_corr_dt, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            count += 1
            if count != 1:
                mt = re.match('(.+)_.+',col[0])
                motifnm = mt.group(1)
                store_motif_gene_line_dic[motifnm] = eachline
            else:
                for i in range(len(col)):
                    geneid = col[i]
                    store_gene_order_dic[geneid] = i + 1

    ##find the most match motif name
    store_final_line_list = []
    store_final_gene_dic = {}
    for eachmotif in store_motif_gene_dic:
        gene_dic = store_motif_gene_dic[eachmotif]
        gene_str = ','.join(list(gene_dic.keys()))

        ##now we will extract information from opt_corr_dt
        if len(list(gene_dic.keys())) != 1:

            gene_line = store_motif_gene_line_dic[eachmotif]
            col = gene_line.strip().split()

            store_gene_val_dic = {}
            gene_val_list = []
            for eachgene in gene_dic:

                if eachgene in store_gene_order_dic:
                    gene_value_order = store_gene_order_dic[eachgene]
                    gene_val = col[gene_value_order]

                    store_gene_val_dic[eachgene] = float(gene_val)
                    gene_val_list.append(gene_val)

                    store_final_gene_dic[eachgene] = 1

                else:
                    print(eachgene + ' not detected for the norm acc')

            if store_gene_val_dic != {}:

                ##updating 122422
                if input_open_compare_corr_rate == 'no':

                    max_gene = max(store_gene_val_dic, key=store_gene_val_dic.get)

                    final_line = eachmotif + '\t' + max_gene + '\t' + str(store_gene_val_dic[max_gene]) + '\t' + gene_str + '\t' + ','.join(gene_val_list)

                else:

                    temp_max_gene = max(store_gene_val_dic, key=store_gene_val_dic.get)

                    ##updating 053123 this is for the 0.1 previous is 0.2
                    if store_gene_val_dic[temp_max_gene] >= 0.1:

                        final_line = eachmotif + '\t' + temp_max_gene + '\t' + str(store_gene_val_dic[temp_max_gene]) + '\t' + gene_str + '\t' + ','.join(gene_val_list)

                    else:

                        store_abs_gene_val_dic = {}
                        for eachabsgene in store_gene_val_dic:
                            store_abs_gene_val_dic[eachabsgene] = abs(store_gene_val_dic[eachabsgene])

                        abs_max_gene = max(store_abs_gene_val_dic, key=store_abs_gene_val_dic.get)

                        final_line = eachmotif + '\t' + abs_max_gene + '\t' + str(store_gene_val_dic[abs_max_gene]) + '\t' + gene_str + '\t' + ','.join(gene_val_list)


                store_final_line_list.append(final_line)



        else:

            gene_line = store_motif_gene_line_dic[eachmotif]
            col = gene_line.strip().split()
            if gene_str in store_gene_order_dic:
                gene_value_order = store_gene_order_dic[gene_str]
                gene_val = col[gene_value_order]

                final_line = eachmotif + '\t' + gene_str + '\t' + gene_val + '\t' + gene_str + '\t' + 'na'
                store_final_line_list.append(final_line)

                store_final_gene_dic[gene_str] = 1
            else:
                print(gene_str + ' not detected for the norm acc')

    with open(input_output_dir + '/opt_corresponding_motif_TF.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ##updating 031323
    ##set all the homo together
    store_final_line_list = []
    with open (input_output_dir + '/opt_corresponding_motif_TF.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            genenmstr = col[3]
            genenm_list = genenmstr.split(',')
            corr_list = col[4].split(',')
            order = -1
            if len(corr_list) >= 2:
                if len(corr_list) == len(genenm_list):
                    for eachgene in genenm_list:
                        order += 1
                        final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + eachgene + '\t' + corr_list[order]
                        store_final_line_list.append(final_line)
            else:
                final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[3] + '\t' + col[2]
                store_final_line_list.append(final_line)

    with open(input_output_dir + '/opt_corresponding_motif_TF_KeepAllGenes.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



    ##filter the opt_corr_dt.txt
    store_final_line_list = []
    count = 0
    with open (input_output_dir + '/opt_corr_dt.txt','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split()
                if col[0] in store_final_gene_dic:
                    store_final_line_list.append(eachline)
            else:
                store_final_line_list.append(eachline)

    with open(input_output_dir + '/opt_corr_dt_flt.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')



corr_TFacc_motifdev (input_blasttf_fl,input_Atgene_TFid_fl,input_Atgene_motif_fl,input_norm_gene_act_fl,input_norm_motif_act_fl,input_output_dir,input_corr_R_script_fl,input_open_compare_corr_rate)





