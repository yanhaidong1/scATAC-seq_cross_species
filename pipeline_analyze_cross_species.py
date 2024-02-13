#!/usr/bin/env python

##
##updating 011724 we will add the methylation to check the LTR
##updating 011624 we will conduct the tf motif enrichment on the specific TE family
##updating 010324 we will update the ss in the step03
##updating 123123 we will add the species-specific ones
##updating 122723 we will add the maize exp for the genes and the common name of rice genes
##updating 122723 we will check how many snps overlapping ACRs for step02
##updating 120723 we will check the cns number per ACR
##udpating 120423 we will add more information about the cell-type-specific
##updating 120423 we will plot the target family motif TF
##updating 113023 we will use the model way to check the motif enrichment
##updating 112723 we will try to check motif enrichment
##updating 112723 we will try to check the conservation for the H3K27me3
##udpating 112123 we will find the case to plot in the browser
##updating 112023 we will add the pm uf to be the blasting
##updating 111823 we will check the motif coverage
##updating 111723 we will update the rice atlas coverage for all syntenic peaks
##updating 111623 we will update df
##updating 111523 we will add the gene cate
##updating 111023 we will add the geneCate to the ACRs
##updating 110923 we will check the overlappings of ACR in all three species
##updating 110823 we will intersect the phylo score to different cate
##updating 110723 we will add if the H3K27me3 overlapping with the rice and maize ACR
##updating 110623 we will also open the sorghum
##updating 110623 we will the different cate of H3K27me3 overlapping ACR
##updating 110323 we will generate different peaks that will be plotted in the browser
##updating 110123 we will add summary of the results

##updating 102623 we will check the syntenic region
##A1 A2 B1 B2
##we will consider the overlapped with the A1 and B2 other than between the A2 and B1

##updating 102023 we will only consider the dic
# ##updating 101723 we will add the cell type specific information to the step01 results
##updating 101523

import re
import glob
import sys
import subprocess
import os
from multiprocessing import Pool
import numpy as np
import os.path
import scipy.stats as stats

from input_required_scripts import s4_subfunctions
from input_required_scripts import s3_subfunctions
from input_required_scripts import s2_subfunctions
from input_required_scripts import s2_subfunctions_exp
from input_required_scripts import s2_subfunctions_cns
from input_required_scripts import s2_subfunctions_snp
from input_required_scripts import s2_subfunctions_nonsyn_syn
from input_required_scripts import s2_subfunctions_H3K27me3
from input_required_scripts import s1_subfunctions
from input_required_scripts import s0_subfunctions
from input_required_scripts import s0_subfunctions_avg_motifs
from input_required_scripts import s0_subfunctions_motifs_TFs
from input_required_scripts import s0_subfunctions_plot_TFs_motifs
from input_required_scripts import s0_subfunctions_motif_model_enrich
from input_required_scripts import s0_subfunctions_cns
from input_required_scripts import s0_subfunctions_nonsyn_syn


##########################
input_required_scripts_dir = sys.argv[1]

########
##step00
input_rice_syntenic_genes_all_os_ACRs_fl = sys.argv[2]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_regions_101723/Os.syntenic_genes.all_os_ACRs.bed

input_all_syntenic_regions_os_acrs_dir = sys.argv[3]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_regions_final_103123/input_all_syntenic_regions_os_acrs_dir


#input_maize_syntenic_regions_os_acrs_fl = sys.argv[3]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_regions_101723/Zm.syntenic_regions.os_acrs.bed

#input_sorghum_syntenic_regions_os_acrs_fl = sys.argv[4]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_regions_101723/Sb.syntenic_regions.os_acrs.bed

##updating 123123
input_rice_acr_add_celltype_fl = sys.argv[4]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_6_characterize_ACRs_060322/output_dir_v8.2_update_051323/step01_ACR_summary/s1_loc_analysis/opt_ACR_sorted_rmdup_addallcelltype.txt


input_total_rice_atlas_acr_fl = sys.argv[5]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_6_characterize_ACRs_060322/output_dir_v8.2_update_051323/step01_ACR_summary/s1_loc_analysis/opt_ACR_sorted_rmdup.txt


##these are the basic functions for the figure 2
input_all_spe_TE_dir = sys.argv[6]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_TE_111323


##this is also used in the step01
##import the rice acr file
##add the cell type specific results
input_all_celltype_acr_dir = sys.argv[7]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_ct_acr_Pablo_110723/final_collect_celltype_acr_dir_111623

##this is also used in the step01
input_all_spe_gene_gff_dir = sys.argv[8]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_gene_gff_111623

input_all_genome_size_dir = sys.argv[9]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_genome_size_111823

input_all_motif_fimo_dir = sys.argv[10]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_motif_fimo_111823

#input_rice_acr_fl = sys.argv[12]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_rice_species_Pablo_091923/Os_data/os_dir/os_acr_classification.no_exons.all_ACRs.classified.sorted.bed

##updating 110123
#input_maize_acr_fl = sys.argv[13]

##updating 110623
#input_sorghum_acr_fl = sys.argv[14]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_maize_sorghum_species_Pablo_091523/sb_acr_classification/sb_acr_classification.no_exons.all_ACRs.classified.sorted.bed

##updating 102623
#input_rice_gff_fl = sys.argv[14]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/Osativa323v7MtPt_sorted.gff3

#input_maize_gff_fl = sys.argv[15]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_maize_sorghum_species_Pablo_091523/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1_noscaf.gff3

#input_sorghum_gff_fl = sys.argv[16]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_maize_sorghum_species_Pablo_091523/Sbicolor.v5.1.prelim.annot/Sbicolorv5.1.gene.gff3


##updating 112723
##check the TFs with motifs
input_gene_raw_sparse_fl_dir = sys.argv[11]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_gene_raw_sparse_112723

input_all_spe_meta_fl_dir = sys.argv[12]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_meta_112723

input_all_spe_SVD_fl_dir = sys.argv[13]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_SVD_112723


input_all_spe_motif_dev_score_dir = sys.argv[14]
#/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_motif_deviation_112723

input_all_spe_ortho_TF_to_AtTF_dir = sys.argv[15]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_to_AtTFs_112723

input_At_TF_to_motif_fl = sys.argv[16]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/11_TFmotif_analysis_040522/output_dir_051523/step02_motif_TF_corresponding_chromVAR/s4_corresponding_TF/s2_correspond_motif_gene_dir/opt_Atgene_motif.txt

input_Atgene_TFid_fl = sys.argv[17]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_cross_species_091823/add_make_At_TF_list_112723/opt_AT_TF_genename_famID_order.txt

input_all_spe_peak_sparse_dir = sys.argv[18]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_022123/add_05_build_PerM_ACR_acc_other_species_091723/01_make_peak_sparse_091723

##updating 120423
ipt_motif_commonnm_fl = sys.argv[19]
##ipt_target_motifs_list_all


########
##step01
input_rice_to_allspe_blast_dir = sys.argv[20]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/collect_rice_to_allspe_blast_dir_112023

input_rice_to_allspe_blast_record_riceID_dir = sys.argv[21]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/collect_rice_to_allspe_blast_record_riceID_dir_112023


#input_rice_to_maize_blast_fl = sys.argv[11]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_blast_final_102623/output_dir/Os_ACRs.vs.Zm.blast_passing_regions.intersecting_regions.bed
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_blast_101623/Os_ACRs.vs.Zm.regions.merged.filtered.blast.blast_passing_regions.intersecting_regions.bed

#input_rice_to_maize_blast_record_riceID_fl = sys.argv[12]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_blast_101623/Os_ACRs.vs.Zm.regions.merged.filtered.blast.blast_passing_regions.intersecting_regions.ref.bed

#input_rice_to_sorghum_blast_fl = sys.argv[13]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_blast_101623/Os_ACRs.vs.Sb.regions.merged.filtered.blast.blast_passing_regions.intersecting_regions.bed

#input_rice_to_sorghum_blast_record_riceID_fl = sys.argv[14]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_blast_101623/Os_ACRs.vs.Sb.regions.merged.filtered.blast.blast_passing_regions.intersecting_regions.ref.bed

input_all_H3K27me3_dir = sys.argv[22]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_022123/collect_all_H3K27me3_fl_dir_112023

#input_rice_H3K27me3_fl = sys.argv[12]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_022123/output_dir_fimo_motif_overlap_riceleaf_InBodyFlank150bpVsNotInH3K27_093023/step01_ACR_H3K27me3_dir/s1_s2_open_add_flank_cate_ACR_H3K27me3_dir/opt_acr_in_and_not_in_H3K27peak_add_ct_addnum_addflank150bp_addcombineCate.txt

#input_maize_H3K27me3_fl = sys.argv[13]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_022123/output_dir_fimo_motif_overlap_maizeleaf_InBodyFlank150bpVsNotInH3K27_093023/step01_ACR_H3K27me3_dir/s1_s2_open_add_flank_cate_ACR_H3K27me3_dir/opt_acr_in_and_not_in_H3K27peak_add_ct_addnum_addflank150bp_addcombineCate.txt

#input_sorghum_H3K27me3_fl = sys.argv[14]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_022123/output_dir_fimo_motif_overlap_sorghumleaf_InBodyFlank150bpVsNotInH3K27_093023/step01_ACR_H3K27me3_dir/s1_s2_open_add_flank_cate_ACR_H3K27me3_dir/opt_acr_in_and_not_in_H3K27peak_add_ct_addnum_addflank150bp_addcombineCate.txt


##updating 110223
#input_rice_acr_celltype_decidedbyCover_fl = sys.argv[15]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_022123/output_dir_fimo_motif_overlap_riceleaf_InBodyFlank150bpVsNotInH3K27_093023/step01_ACR_H3K27me3_dir/s1_s2_open_add_flank_cate_ACR_H3K27me3_dir/opt_acr_in_and_not_in_H3K27peak_add_ct_addnum_addflank150bp_addcombineCate.txt

#input_maize_acr_celltype_decidedbyCover_fl = sys.argv[16]

##updating 110623
#input_sorghum_acr_celltype_decidedbyCover_fl = sys.argv[17]
##


input_all_syntenic_genes_all_os_ACRs_dir = sys.argv[23]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_regions_final_103123/input_all_syntenic_genes.all_os_ACRs_dir


#input_maize_syntenic_genes_all_os_ACRs_fl = sys.argv[15]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_regions_final_102023/Zm.syntenic_genes.all_os_ACRs.bed

#input_sorghum_syntenic_genes_all_os_ACRs_fl = sys.argv[16]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_five_species_Pablo_syntenic_block_092223/syntenicBlocks_regions_final_103123/Sb.syntenic_genes.all_os_ACRs.bed






##updating 110123
########
##step02
input_rice_phy_score_fl = sys.argv[24]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_07_check_conservation_score_052923/output_dir/step02_analyze_conservation_scores_dir/temp_conservation_score_sorted.txt
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_conservation_score_phyloP_dir/Osativa_phyloP.renamed.bdg

##updating 111723
input_rice_atlas_acr_coverage_fl = sys.argv[25]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_5_iden_DApeaks_040422/02_pipeline_callDAR_040522/output_dir_CoverageVer_newCalling_051323/step02_add_peaks_overlap_positive_negative_perorganct_dir/temp_organcelltype_peak_largest_readCover_middle_addPosiNeg.txt


input_all_spe_gene_annot_fl_dir = sys.argv[26]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_gene_annot_dir_120623


##updating 120723
input_all_spe_cns_gff_fl_dir = sys.argv[27]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_cns_gff_120723

##for the Pablo version
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_CNS_Pablo_121823/Pablo_CNS_all_species_to_others_gff_dir


input_marker_gene_list_fl_dir = sys.argv[28]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_marker_gene_120723

ipt_cpm_snRNAseq_fl = sys.argv[29]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/05_integrate_snRNAseq_041322/01_Data_process_090922/04_annotate_by_scATACseq_data_042623/output_dir/step02_annotate_by_correlation_acc_exp_dir/opt_seedlingsnRNAseq_perM_genes_exp_clusters.txt

ipt_cluster_annot_fl = sys.argv[30]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_022123/input_seedling_snRNAseq_markers_dir/ipt_cluster_annotation_fl_seedling.txt

#3updating 122723
ipt_maize_cpm_snRNAseq_fl = sys.argv[31]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_022123/add_04_build_PerM_gene_expression_other_species_091723/output_dir/opt_maizeLeaf_perM_genes_exp_clusters.txt

input_all_spe_geneID_symbol_dir = sys.argv[32]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_gene_symbol_122723

input_all_spe_ortho_dir = sys.argv[33]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_five_species_ortho_genes_122823

##updating 122723
input_rice_3k_snp_fl = sys.argv[34]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_rice_snp_122723/NB_bialSNP_pseudo_canonical_ALL.vcf


##updating 010524
input_H3K27me3_fl = sys.argv[35]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_022123/input_H3K27_bw_fl_dir/opt_modi_chr_bw_dir/Eseedling_H3K27me3.bw


##updating 011724
input_DNAmethylation_bw_dir = sys.argv[36]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_7_methylation_ACRs_062422/collect_bw_fl_022223/leafmethy


########
##step03

input_target_organ_H3K27_peak_all_dir = sys.argv[37]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_03_more_chipSeq_data_030223/05_make_final_peaks_for_diff_species_091523/collect_all_H3K27me3_final_peaks_112023


#ipt_target_organ_H3K27_peak_rice_fl = sys.argv[18]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_03_more_chipSeq_data_030223/05_make_final_peaks_for_diff_species_091523/output_dir_rice_leaf_H3K27me3_peaks_092023/store_epic2_fixby_macs2_dir/opt_RiceLeafH3K27me3_final_peak_sorted.txt

#ipt_target_organ_H3K27_peak_maize_fl = sys.argv[19]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_03_more_chipSeq_data_030223/05_make_final_peaks_for_diff_species_091523/output_dir_maize_leaf_H3K27me3_peaks_091623/store_epic2_fixby_macs2_dir/opt_MaizeLeafH3K27me3_final_peak_sorted.txt

##updating 110623
#ipt_target_organ_H3K27_peak_sorghum_fl = sys.argv[20]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_03_more_chipSeq_data_030223/05_make_final_peaks_for_diff_species_091523/output_dir_sorghum_leaf_H3K27me3_peaks_091823/store_epic2_fixby_macs2_dir/opt_SorghumLeafH3K27me3_final_peak_sorted.txt

input_rice_all_organ_H3K27_fl_dir = sys.argv[38]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/store_all_organ_rice_H3K27me3_dir_112823/

##updating 121223
input_all_spe_PRE_motif_fimo_dir = sys.argv[39]
##/scratch/hy17471/rice_altlas_scATAC_seq_042021/07_4_callpeaks_FDRway_pipelineVer_040322/add_02_check_H3K27_motifs_022123/collect_all_motif_fimo_PRE_pval00002_dir_121223

input_configure_fl = sys.argv[40]

input_output_dir = sys.argv[41]





########
##step00
##we will make a summary of how many syntenic regions

##1 How many rice syntenic regions, and what's the size of the syntenic block to the whole genome
##2 what's the distribution of distance per syntenic regions
##3 How many rice ACR in the syntenic regions (Here we would find the not shared ACRs in the rice)
##4 add these not shared ACRs to the step01
##5 step01 need to add the broad or restricted ones


def step00_basic_characters (input_required_scripts_dir,input_output_dir,input_rice_syntenic_genes_all_os_ACRs_fl,
                             input_all_syntenic_regions_os_acrs_dir,
                             input_all_H3K27me3_dir,input_rice_acr_add_celltype_fl,input_total_rice_atlas_acr_fl,
                             input_all_spe_TE_dir,input_all_celltype_acr_dir,input_all_spe_gene_gff_dir,input_all_genome_size_dir,
                             input_all_motif_fimo_dir,
                             input_gene_raw_sparse_fl_dir,input_all_spe_meta_fl_dir,input_all_spe_SVD_fl_dir,
                             input_all_spe_motif_dev_score_dir,
                             input_all_spe_ortho_TF_to_AtTF_dir,input_At_TF_to_motif_fl,input_Atgene_TFid_fl,
                             input_all_spe_peak_sparse_dir,ipt_motif_commonnm_fl,
                             input_all_spe_cns_gff_fl_dir,
                             input_core_num,
                             s0_s1_open_check_syntenic_block,s0_s1_target_spe1,s0_s1_target_spe2_str,s0_s1_width_bin,
                             s0_s2_open_check_te_composition,
                             s0_s3_open_check_gene_cate,
                             s0_s4_open_motif_coverage,s0_s4_target_species_str,s0_s4_sub_open_build_control,
                             s0_s4_sub_open_check_avg_for_acr,s0_s4_sub_extend_acr_range_bp,
                             s0_s4_sub_open_build_posi_coverage,
                             s0_s5_open_motifs_TFs,
                             s0_s5_sub_open_prepare_gene_accessibility,s0_s5_sub_target_spe_str_for_geneacc_smooth,
                             s0_s5_sub_open_smooth_gene,s0_s5_sub_target_cluster,
                             s0_s5_sub_open_smooth_motif,
                             s0_s5_sub_open_corr_motif_TF,
                             s0_s5_sub_open_plot_TF_motif,s0_s5_sub_lim, s0_s5_sub_plot_topnum_pairs,
                             s0_s5_sub_open_plot_TFfam_TF_motif,s0_s5_sub_target_fam_str,
                             s0_s5_sub_open_plot_high_corr_all_spe,s0_s5_sub_target_plot_spe_str, s0_s5_sub_target_plot_corr_cutoff,
                             s0_s5_sub_open_plot_Uf_TF_motif,s0_s5_sub_motif_TF_str_for_Uf,
                             s0_s6_model_based_motif_enrich,s0_s6_sub_open_flt_meta_cell,s0_s6_sub_cluster_colnum_str,
                             s0_s6_sub_target_spe_str_for_add_ACRnum,
                             s0_s6_sub_open_add_ACRnum_to_meta,
                             s0_s6_sub_open_add_log10tn5_lib_to_meta,s0_s6_sub_target_library_colnum,s0_s6_sub_target_tn5_colnum,
                             s0_s6_sub_open_prepare_acrmotif,s0_s6_sub_target_spe_str_for_conduct_enrich_str,
                             s0_s6_sub_open_prepare_enrich,s0_s6_sub_pvalcutoff,
                             s0_s6_sub_open_conduct_enrich,s0_s6_sub_target_clusternm,
                             s0_s7_open_check_cns,s0_s7_sub_open_eachspe_overlapCNS,s0_s7_sub_target_all_spe_str,s0_s7_sub_open_rice_atlas_CNS):

    step00_basic_characters_dir = input_output_dir + '/step00_basic_characters_dir'
    if not os.path.exists(step00_basic_characters_dir):
        os.makedirs(step00_basic_characters_dir)

    if s0_s1_open_check_syntenic_block == 'yes':

        s0_s1_open_check_syntenic_block_dir = step00_basic_characters_dir + '/s0_s1_open_check_syntenic_block_dir'
        if not os.path.exists(s0_s1_open_check_syntenic_block_dir):
            os.makedirs(s0_s1_open_check_syntenic_block_dir)

        s0_s1_target_spe2_list = s0_s1_target_spe2_str.split(',')

        ##updating 112023
        ##we will use the loop other than run several times
        for eachspe2 in s0_s1_target_spe2_list:


            ##for the rice to maize
            store_rice_to_maize_dir = s0_s1_open_check_syntenic_block_dir + '/store_' + s0_s1_target_spe1 + '_to_' + eachspe2 + '_dir'
            if not os.path.exists(store_rice_to_maize_dir):
                os.makedirs(store_rice_to_maize_dir)


            spe1_prefix = s0_s1_target_spe1
            spe2_prefix = eachspe2

            input_maize_syntenic_regions_os_acrs_fl = input_all_syntenic_regions_os_acrs_dir + '/' + spe2_prefix + '.syntenic_regions.os_acrs.bed'

            input_rice_H3K27me3_fl = input_all_H3K27me3_dir + '/' + spe1_prefix + '_H3K27me3_cate.txt'

            s0_subfunctions.subfunction_check_syntenic_blocks(input_rice_syntenic_genes_all_os_ACRs_fl, input_maize_syntenic_regions_os_acrs_fl, input_rice_H3K27me3_fl,
                                              store_rice_to_maize_dir,
                                              spe1_prefix, spe2_prefix)

        ##for the rice to build the all the syntenic region file
        spe1_prefix = s0_s1_target_spe1
        s0_subfunctions.subfunction_build_syntenic_bed(input_rice_syntenic_genes_all_os_ACRs_fl, input_total_rice_atlas_acr_fl,s0_s1_open_check_syntenic_block_dir, spe1_prefix)


        ##updating 123123
        ##we will check the syntenic region for atlas ACR
        input_rice_atlas_acr_fl = input_total_rice_atlas_acr_fl
        ##s0_s1_open_check_syntenic_block_dir/opt_allrice_syntenic_region_add_distance_sorted.txt
        ipt_syntenic_block_fl = s0_s1_open_check_syntenic_block_dir + '/opt_allrice_syntenic_region_add_distance_sorted.txt'
        s0_subfunctions_nonsyn_syn.subfunction_check_overlap_number_ACRs_in_syn_non_syn(input_rice_acr_add_celltype_fl, input_rice_atlas_acr_fl,
                                                             ipt_syntenic_block_fl, s0_s1_open_check_syntenic_block_dir, s0_s1_width_bin)

        s0_subfunctions_nonsyn_syn.subfunction_check_size_of_syntenic_region(s0_s1_open_check_syntenic_block_dir)



    if s0_s2_open_check_te_composition == 'yes':

        s0_s2_open_check_te_composition_dir = step00_basic_characters_dir + '/s0_s2_open_check_te_composition_dir'
        if not os.path.exists(s0_s2_open_check_te_composition_dir):
            os.makedirs(s0_s2_open_check_te_composition_dir)

        s0_s2_target_species_list = s0_s1_target_spe2_str.split(',')

        all_fl_list = glob.glob(input_all_spe_TE_dir + '/*')
        for eachfl in all_fl_list:
            mt = re.match('.+/(.+)',eachfl)
            flnm = mt.group(1)

            mt = re.match('(.+)\..+\.fa\.mod\.EDTA\.TEanno\.gff3',flnm)
            spnm = mt.group(1)

            if spnm in s0_s2_target_species_list:

                spnm_dir = s0_s2_open_check_te_composition_dir + '/' + spnm
                if not os.path.exists(spnm_dir):
                    os.makedirs(spnm_dir)

                ipt_acr_fl = input_all_celltype_acr_dir + '/'+ spnm + '_acr_celltype.txt'

                s0_subfunctions.subfunction_check_TE_composition(eachfl, ipt_acr_fl, spnm, spnm_dir)


    if s0_s3_open_check_gene_cate == 'yes':

        s0_s3_open_check_gene_cate_dir = step00_basic_characters_dir + '/s0_s3_open_check_gene_cate_dir'
        if not os.path.exists(s0_s3_open_check_gene_cate_dir):
            os.makedirs(s0_s3_open_check_gene_cate_dir)

        all_spe_fl_list = glob.glob(input_all_spe_gene_gff_dir + '/*')
        for eachfl in all_spe_fl_list:

            mt = re.match('.+/(.+)',eachfl)
            flnm = mt.group(1)
            mt = re.match('(.+)_gene\.gff',flnm)
            spnm = mt.group(1)

            opt_spm_dir = s0_s3_open_check_gene_cate_dir + '/' + spnm
            if not os.path.exists(opt_spm_dir):
                os.makedirs(opt_spm_dir)


            input_gff_fl = input_all_spe_gene_gff_dir + '/' + spnm + '_gene.gff'
            input_acr_fl = input_all_celltype_acr_dir + '/' + spnm + '_acr_celltype.txt'

            s0_subfunctions.subfunction_make_acr_close_to_gene_cate(input_gff_fl, input_acr_fl,
                                                                    opt_spm_dir, spnm)


    if s0_s4_open_motif_coverage == 'yes':

        s0_s4_open_motif_coverage_dir = step00_basic_characters_dir + '/s0_s4_open_motif_coverage_dir'
        if not os.path.exists(s0_s4_open_motif_coverage_dir):
            os.makedirs(s0_s4_open_motif_coverage_dir)


        if s0_s4_sub_open_build_control == 'yes':

            s0_s4_sub_open_build_control_dir = s0_s4_open_motif_coverage_dir + '/s0_s4_sub_open_build_control_dir'
            if not os.path.exists(s0_s4_sub_open_build_control_dir):
                os.makedirs(s0_s4_sub_open_build_control_dir)

            s0_s4_target_species_list = s0_s4_target_species_str.split(',')

            for eachspe in s0_s4_target_species_list:

                opt_spm_dir = s0_s4_sub_open_build_control_dir + '/' + eachspe
                if not os.path.exists(opt_spm_dir):
                    os.makedirs(opt_spm_dir)

                ipt_spe_te_gff_fl = ''
                all_fl_list = glob.glob(input_all_spe_TE_dir + '/*')
                for eachfl in all_fl_list:
                    mt = re.match('.+/(.+)', eachfl)
                    flnm = mt.group(1)
                    mt = re.match('(.+)\..+\.fa\.mod\.EDTA\.TEanno\.gff3', flnm)
                    spnm = mt.group(1)
                    if spnm == eachspe:
                        ipt_spe_te_gff_fl = eachfl

                if ipt_spe_te_gff_fl != '':
                    print('the te file is ' + ipt_spe_te_gff_fl)
                else:
                    print('Please check if the TE file has been identified')
                    break

                ipt_spe_acr_fl = input_all_celltype_acr_dir + '/'+ eachspe + '_acr_celltype.txt'
                ipt_spe_genome_size_fl = input_all_genome_size_dir + '/' + eachspe + '_genome_size.txt'


                s0_subfunctions_avg_motifs.subfunction_s1_create_control(ipt_spe_te_gff_fl, ipt_spe_acr_fl, ipt_spe_genome_size_fl, opt_spm_dir)


        if s0_s4_sub_open_check_avg_for_acr == 'yes':

            s0_s4_sub_open_check_avg_for_acr_dir = s0_s4_open_motif_coverage_dir + '/s0_s4_sub_open_check_avg_for_acr_dir'
            if not os.path.exists(s0_s4_sub_open_check_avg_for_acr_dir):
                os.makedirs(s0_s4_sub_open_check_avg_for_acr_dir)

            s0_s4_target_species_list = s0_s4_target_species_str.split(',')

            for eachspe in s0_s4_target_species_list:

                opt_spm_dir = s0_s4_sub_open_check_avg_for_acr_dir + '/' + eachspe
                if not os.path.exists(opt_spm_dir):
                    os.makedirs(opt_spm_dir)

                true_dir = opt_spm_dir + '/' + 'true_dir'
                if not os.path.exists(true_dir):
                    os.makedirs(true_dir)

                ipt_spe_acr_fl = input_all_celltype_acr_dir + '/' + eachspe + '_acr_celltype.txt'

                ipt_spe_motif_fl = input_all_motif_fimo_dir + '/' + eachspe + '_motif_fimo.txt'

                ipt_convert2matrix_script = input_required_scripts_dir + '/convert2matrix.pl'
                s0_subfunctions_avg_motifs.subfunction_s2_check_avg_for_true_acr(true_dir, ipt_spe_acr_fl, ipt_spe_motif_fl, ipt_convert2matrix_script,
                                                                                s0_s4_sub_extend_acr_range_bp)

                control_dir = opt_spm_dir + '/' + 'control_dir'
                if not os.path.exists(control_dir):
                    os.makedirs(control_dir)

                ipt_control_acr_fl = s0_s4_open_motif_coverage_dir + '/s0_s4_sub_open_build_control_dir' + '/' + eachspe + '/opt_candidate_control_regions.txt'


                s0_subfunctions_avg_motifs.subfunction_s2_check_avg_for_control_acr (control_dir,ipt_control_acr_fl,ipt_spe_motif_fl,s0_s4_sub_extend_acr_range_bp,
                                           ipt_convert2matrix_script)

        if s0_s4_sub_open_build_posi_coverage == 'yes':

            s0_s4_sub_open_build_posi_coverage_dir = s0_s4_open_motif_coverage_dir + '/s0_s4_sub_open_build_posi_coverage_dir'
            if not os.path.exists(s0_s4_sub_open_build_posi_coverage_dir):
                os.makedirs(s0_s4_sub_open_build_posi_coverage_dir)

            s0_s4_target_species_list = s0_s4_target_species_str.split(',')

            for eachspe in s0_s4_target_species_list:

                opt_spm_dir = s0_s4_sub_open_build_posi_coverage_dir + '/' + eachspe
                if not os.path.exists(opt_spm_dir):
                    os.makedirs(opt_spm_dir)

                ipt_position_count_mtx_true_acr_fl = s0_s4_open_motif_coverage_dir + '/s0_s4_sub_open_check_avg_for_acr_dir/' + \
                                                     eachspe + '/true_dir/opt_position_count_mtx_true_acr.txt'
                ipt_position_count_mtx_control_acr_fl = s0_s4_open_motif_coverage_dir + '/s0_s4_sub_open_check_avg_for_acr_dir/' + \
                                                        eachspe + '/control_dir/opt_position_count_mtx_control_acr.txt'
                ipt_R_script_for_build_coverage = input_required_scripts_dir + '/build_coverage_data.R'
                s0_subfunctions_avg_motifs.subfunction_s3_build_posi_coverage(opt_spm_dir, ipt_position_count_mtx_true_acr_fl,
                                                   ipt_position_count_mtx_control_acr_fl,
                                                   ipt_R_script_for_build_coverage)


    if s0_s5_open_motifs_TFs == 'yes':

        s0_s5_open_motifs_TFs_dir = step00_basic_characters_dir + '/s0_s5_open_motifs_TFs_dir'
        if not os.path.exists(s0_s5_open_motifs_TFs_dir):
            os.makedirs(s0_s5_open_motifs_TFs_dir)

        if s0_s5_sub_open_prepare_gene_accessibility == 'yes':

            ##we will smooth smooth the gene accessibility
            store_prepare_gene_acc_dir = s0_s5_open_motifs_TFs_dir + '/store_prepare_gene_acc_dir'
            if not os.path.exists(store_prepare_gene_acc_dir):
                os.makedirs(store_prepare_gene_acc_dir)

            s0_s5_sub_target_spe_str_for_geneacc_list = s0_s5_sub_target_spe_str_for_geneacc_smooth.split(',')

            for eachspe in s0_s5_sub_target_spe_str_for_geneacc_list:

                store_prepare_gene_acc_eachspe_dir = store_prepare_gene_acc_dir + '/' + eachspe
                if not os.path.exists(store_prepare_gene_acc_eachspe_dir):
                    os.makedirs(store_prepare_gene_acc_eachspe_dir)

                input_gene_cell_raw_sparse_fl = input_gene_raw_sparse_fl_dir + '/' + eachspe + '_raw_gene_sparse.txt'
                input_meta_fl = input_all_spe_meta_fl_dir + '/' + eachspe + '_meta.txt'

                s0_subfunctions_motifs_TFs.subfunction_prepare_gene_accessibility(input_gene_cell_raw_sparse_fl, input_meta_fl,
                                                                                  input_required_scripts_dir, eachspe,
                                                                                  store_prepare_gene_acc_eachspe_dir)

        if s0_s5_sub_open_smooth_gene == 'yes':

            store_smooth_gene_dir = s0_s5_open_motifs_TFs_dir + '/store_smooth_gene_dir'
            if not os.path.exists(store_smooth_gene_dir):
                os.makedirs(store_smooth_gene_dir)

            s0_s5_sub_target_spe_str_for_smooth_list = s0_s5_sub_target_spe_str_for_geneacc_smooth.split(',')

            for eachspe in s0_s5_sub_target_spe_str_for_smooth_list:

                store_smooth_eachspe_dir = store_smooth_gene_dir + '/' + eachspe
                if not os.path.exists(store_smooth_eachspe_dir):
                    os.makedirs(store_smooth_eachspe_dir)

                input_meta_fl = input_all_spe_meta_fl_dir + '/' + eachspe + '_meta.txt'
                input_gene_norm_sparse_fl = s0_s5_open_motifs_TFs_dir + '/store_prepare_gene_acc_dir/' + eachspe + '/genebody_accessibility_dir/' + eachspe + '.GBaccessibility.sparse'
                input_svd_fl = input_all_spe_SVD_fl_dir + '/' + eachspe + '_svd.txt'

                s0_subfunctions_motifs_TFs.subfunction_prepare_smooth_gene_accessibility (input_required_scripts_dir,input_meta_fl,input_gene_norm_sparse_fl,input_svd_fl,
                                                    input_core_num,s0_s5_sub_target_cluster,store_smooth_eachspe_dir)

        if s0_s5_sub_open_smooth_motif == 'yes':

            store_smooth_motif_dir = s0_s5_open_motifs_TFs_dir + '/store_smooth_motif_dir'
            if not os.path.exists(store_smooth_motif_dir):
                os.makedirs(store_smooth_motif_dir)

            s0_s5_sub_target_spe_str_for_smooth_list = s0_s5_sub_target_spe_str_for_geneacc_smooth.split(',')

            for eachspe in s0_s5_sub_target_spe_str_for_smooth_list:

                store_smooth_eachspe_dir = store_smooth_motif_dir + '/' + eachspe
                if not os.path.exists(store_smooth_eachspe_dir):
                    os.makedirs(store_smooth_eachspe_dir)

                input_meta_fl = input_all_spe_meta_fl_dir + '/' + eachspe + '_meta.txt'
                ipt_motif_deviation_fl = input_all_spe_motif_dev_score_dir + '/' + eachspe + '.motif.deviation.txt'
                ipt_svd_fl = input_all_spe_SVD_fl_dir + '/' + eachspe + '_svd.txt'

                s0_subfunctions_motifs_TFs.subfunction_prepare_smooth_motif_deviation(input_required_scripts_dir, input_meta_fl,
                                                           ipt_motif_deviation_fl,
                                                           ipt_svd_fl, eachspe, s0_s5_sub_target_cluster, store_smooth_eachspe_dir)

        if s0_s5_sub_open_corr_motif_TF == 'yes':

            store_corr_motif_TF_dir = s0_s5_open_motifs_TFs_dir + '/store_corr_motif_TF_dir'
            if not os.path.exists(store_corr_motif_TF_dir):
                os.makedirs(store_corr_motif_TF_dir)

            s0_s5_sub_target_spe_str_for_smooth_list = s0_s5_sub_target_spe_str_for_geneacc_smooth.split(',')

            for eachspe in s0_s5_sub_target_spe_str_for_smooth_list:

                store_motif_TF_dir = store_corr_motif_TF_dir + '/' + eachspe
                if not os.path.exists(store_motif_TF_dir):
                    os.makedirs(store_motif_TF_dir)

                ##ipt_smooth_sparse_fl = opt_s3_smooth_motif + '/opt_motif_smooth_deviations.sparse'
                ##ipt_smooth_TF_mtx_rds_fl = step02_s4_3_sub_use_smooth_TF_acc_fl_path

                ipt_motif_smooth_sparse_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_motif_dir/' + eachspe + '/opt_motif_smooth_deviations.sparse'
                ipt_smooth_TF_mtx_rds_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_gene_dir/' + eachspe + '/opt_allgenes_impute.activity.rds'
                input_meta_fl = input_all_spe_meta_fl_dir + '/' + eachspe + '_meta.txt'
                ipt_correspond_At_gene_to_spe_fl = input_all_spe_ortho_TF_to_AtTF_dir + '/' + eachspe + '.orthos_from.At.markers.tis_leaf.visualize.txt'


                s0_subfunctions_motifs_TFs.subfunction_correlate_TF_motif(input_required_scripts_dir, ipt_motif_smooth_sparse_fl, input_meta_fl,
                                                                           ipt_smooth_TF_mtx_rds_fl, ipt_correspond_At_gene_to_spe_fl,
                                                                           input_At_TF_to_motif_fl, input_Atgene_TFid_fl,
                                                                           s0_s5_sub_target_cluster, store_motif_TF_dir)

        if s0_s5_sub_open_plot_TF_motif == 'yes':

            store_plot_TF_motif_dir = s0_s5_open_motifs_TFs_dir + '/store_plot_TF_motif_dir'
            if not os.path.exists(store_plot_TF_motif_dir):
                os.makedirs(store_plot_TF_motif_dir)

            s0_s5_sub_target_spe_str_for_smooth_list = s0_s5_sub_target_spe_str_for_geneacc_smooth.split(',')

            for eachspe in s0_s5_sub_target_spe_str_for_smooth_list:

                store_motif_TF_dir = store_plot_TF_motif_dir + '/' + eachspe
                if not os.path.exists(store_motif_TF_dir):
                    os.makedirs(store_motif_TF_dir)

                ipt_TF_motif_corr_fl = s0_s5_open_motifs_TFs_dir + '/store_corr_motif_TF_dir/' + eachspe + '/opt_corresponding_motif_TF.txt'
                ipt_gene_impute_acc_rds_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_gene_dir/' + eachspe + '/opt_allgenes_impute.activity.rds'
                GAobj_rds_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_gene_dir/' + eachspe + '/GAobj.rds'

                ##plot the TFs
                s0_subfunctions_plot_TFs_motifs.subfunction_plot_TFs(input_required_scripts_dir, ipt_TF_motif_corr_fl, ipt_gene_impute_acc_rds_fl,
                                     GAobj_rds_fl
                                     , store_motif_TF_dir, s0_s5_sub_lim, s0_s5_sub_plot_topnum_pairs)

                ipt_motif_impute_acc_threesparse_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_motif_dir/' + eachspe + '/opt_motif_smooth_deviations.sparse'
                ipt_meta_fl = input_all_spe_meta_fl_dir + '/' + eachspe + '_meta.txt'

                ##plot the motifs
                s0_subfunctions_plot_TFs_motifs.subfunction_plot_motifs(input_required_scripts_dir, ipt_TF_motif_corr_fl, ipt_motif_impute_acc_threesparse_fl,
                                        ipt_meta_fl,store_motif_TF_dir, s0_s5_sub_lim, s0_s5_sub_plot_topnum_pairs)


        ##updating 120423
        if s0_s5_sub_open_plot_TFfam_TF_motif == 'yes':

            store_plot_TFfam_TF_motif_dir = s0_s5_open_motifs_TFs_dir + '/store_plot_TFfam_TF_motif_dir'
            if not os.path.exists(store_plot_TFfam_TF_motif_dir):
                os.makedirs(store_plot_TFfam_TF_motif_dir)

            s0_s5_sub_target_spe_str_for_smooth_list = s0_s5_sub_target_spe_str_for_geneacc_smooth.split(',')

            for eachspe in s0_s5_sub_target_spe_str_for_smooth_list:

                store_motif_TF_dir = store_plot_TFfam_TF_motif_dir + '/' + eachspe
                if not os.path.exists(store_motif_TF_dir):
                    os.makedirs(store_motif_TF_dir)

                ipt_TF_motif_corr_fl = s0_s5_open_motifs_TFs_dir + '/store_corr_motif_TF_dir/' + eachspe + '/opt_corresponding_motif_TF.txt'
                ipt_gene_impute_acc_rds_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_gene_dir/' + eachspe + '/opt_allgenes_impute.activity.rds'
                GAobj_rds_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_gene_dir/' + eachspe + '/GAobj.rds'

                ipt_motif_impute_acc_threesparse_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_motif_dir/' + eachspe + '/opt_motif_smooth_deviations.sparse'
                ipt_meta_fl = input_all_spe_meta_fl_dir + '/' + eachspe + '_meta.txt'

                s0_subfunctions_plot_TFs_motifs.subfunction_plot_TFfam_motifs_TFs(input_required_scripts_dir, ipt_TF_motif_corr_fl,
                                                                                  ipt_motif_commonnm_fl,
                                                                                  ipt_gene_impute_acc_rds_fl, GAobj_rds_fl,
                                                                                  ipt_motif_impute_acc_threesparse_fl, ipt_meta_fl,
                                                                                  store_motif_TF_dir, s0_s5_sub_target_fam_str, s0_s5_sub_lim)


        ##updating 123023
        if s0_s5_sub_open_plot_Uf_TF_motif == 'yes':

            s0_s5_sub_open_plot_Uf_TF_motif_dir = s0_s5_open_motifs_TFs_dir + '/s0_s5_sub_open_plot_Uf_TF_motif_dir'
            if not os.path.exists(s0_s5_sub_open_plot_Uf_TF_motif_dir):
                os.makedirs(s0_s5_sub_open_plot_Uf_TF_motif_dir)

            eachspe = 'Uf'

            store_motif_TF_dir = s0_s5_sub_open_plot_Uf_TF_motif_dir + '/' + 'Uf'
            if not os.path.exists(store_motif_TF_dir):
                os.makedirs(store_motif_TF_dir)


            ipt_gene_impute_acc_rds_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_gene_dir/' + eachspe + '/opt_allgenes_impute.activity.rds'
            GAobj_rds_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_gene_dir/' + eachspe + '/GAobj.rds'

            ipt_motif_impute_acc_threesparse_fl = s0_s5_open_motifs_TFs_dir + '/store_smooth_motif_dir/' + eachspe + '/opt_motif_smooth_deviations.sparse'
            ipt_meta_fl = input_all_spe_meta_fl_dir + '/' + eachspe + '_meta.txt'

            ipt_motif_TF_str = s0_s5_sub_motif_TF_str_for_Uf

            s0_subfunctions_plot_TFs_motifs.subfunction_plot_motifs_for_Uf (input_required_scripts_dir,ipt_motif_TF_str,ipt_motif_impute_acc_threesparse_fl, ipt_meta_fl,
                                                                            ipt_gene_impute_acc_rds_fl,GAobj_rds_fl,
                                                                            store_motif_TF_dir,s0_s5_sub_lim)



        ##updating 121023
        if s0_s5_sub_open_plot_high_corr_all_spe == 'yes':

            s0_s5_sub_open_plot_high_corr_all_spe_dir = s0_s5_open_motifs_TFs_dir + '/s0_s5_sub_open_plot_high_corr_all_spe_dir'
            if not os.path.exists(s0_s5_sub_open_plot_high_corr_all_spe_dir):
                os.makedirs(s0_s5_sub_open_plot_high_corr_all_spe_dir)

            ipt_correlation_dir = s0_s5_open_motifs_TFs_dir + '/store_corr_motif_TF_dir'


            s0_subfunctions_plot_TFs_motifs.subfunction_select_target_TF_motif_in_all_spe_to_plot(ipt_correlation_dir, ipt_motif_commonnm_fl,
                                                                                            s0_s5_sub_target_plot_spe_str, s0_s5_sub_target_plot_corr_cutoff,
                                                                                                  s0_s5_sub_open_plot_high_corr_all_spe_dir)




    if s0_s6_model_based_motif_enrich == 'yes':

        s0_s6_model_based_motif_enrich_dir = step00_basic_characters_dir + '/s0_s6_model_based_motif_enrich_dir'
        if not os.path.exists(s0_s6_model_based_motif_enrich_dir):
            os.makedirs(s0_s6_model_based_motif_enrich_dir)
        

        #ipt_spe_motif_fl = input_all_motif_fimo_dir + '/' + eachspe + '_motif_fimo.txt'

        if s0_s6_sub_open_flt_meta_cell == 'yes':

            store_filter_cell_meta_dir = s0_s6_model_based_motif_enrich_dir + '/store_filter_cell_meta_dir'
            if not os.path.exists(store_filter_cell_meta_dir):
                os.makedirs(store_filter_cell_meta_dir)


            ##we will first check the meta data
            s0_subfunctions_motif_model_enrich.subfunction_flt_meta_cells_keep_same_cellnum_celltypes (input_all_spe_meta_fl_dir,
                                                                                                       s0_s6_sub_cluster_colnum_str,
                                                                                                       store_filter_cell_meta_dir)

        if s0_s6_sub_open_add_ACRnum_to_meta == 'yes':

            store_add_ACRnum_to_meta_dir = s0_s6_model_based_motif_enrich_dir + '/store_add_ACRnum_to_meta_dir'
            if not os.path.exists(store_add_ACRnum_to_meta_dir):
                os.makedirs(store_add_ACRnum_to_meta_dir)

            s0_s6_sub_target_spe_str_list = s0_s6_sub_target_spe_str_for_add_ACRnum.split(',')

            for eachspe in s0_s6_sub_target_spe_str_list:

                ipt_meta_fl = s0_s6_model_based_motif_enrich_dir + '/store_filter_cell_meta_dir' + '/' + eachspe + '_flt_meta.txt'
                ipt_peak_sparse_fl = input_all_spe_peak_sparse_dir + '/' + 'output_dir_' + eachspe \
                                     + '/opt_peak_sparse_dir/opt_peak_' + eachspe + 'Leaf.sparse'

                s0_subfunctions_motif_model_enrich.subfunction_build_peak_cell_in_meta(ipt_meta_fl,
                                                                                       ipt_peak_sparse_fl,
                                                                                       store_add_ACRnum_to_meta_dir, eachspe)

        if s0_s6_sub_open_add_log10tn5_lib_to_meta == 'yes':

            store_add_log10tn5_lib_dir = s0_s6_model_based_motif_enrich_dir + '/store_add_log10tn5_lib_dir'
            if not os.path.exists(store_add_log10tn5_lib_dir):
                os.makedirs(store_add_log10tn5_lib_dir)

            s0_s6_sub_target_spe_str_list = s0_s6_sub_target_spe_str_for_add_ACRnum.split(',')

            for eachspe in s0_s6_sub_target_spe_str_list:

                ipt_meta_fl = s0_s6_model_based_motif_enrich_dir + '/store_add_ACRnum_to_meta_dir/' + '/opt_' + eachspe + '.txt'

                s0_subfunctions_motif_model_enrich.subfunction_add_log10tn5_library(ipt_meta_fl, s0_s6_sub_target_library_colnum,
                                                                                    s0_s6_sub_target_tn5_colnum, eachspe,
                                                                                    store_add_log10tn5_lib_dir)

        if s0_s6_sub_open_prepare_acrmotif == 'yes':

            store_prepare_arcmotif_dir = s0_s6_model_based_motif_enrich_dir + '/store_prepare_arcmotif_dir'
            if not os.path.exists(store_prepare_arcmotif_dir):
                os.makedirs(store_prepare_arcmotif_dir)

            s0_s6_sub_target_spe_str_for_conduct_enrich_list = s0_s6_sub_target_spe_str_for_conduct_enrich_str.split(',')

            for eachspe in s0_s6_sub_target_spe_str_for_conduct_enrich_list:

                ipt_spe_ACR_fl =  input_all_celltype_acr_dir + '/' + eachspe + '_acr_celltype.txt'

                ipt_motif_fl = input_all_motif_fimo_dir + '/' + eachspe + '_motif_fimo.txt'

                s0_subfunctions_motif_model_enrich.subfunction_intersect_acr_motif(ipt_spe_ACR_fl, ipt_motif_fl,
                                                                                   store_prepare_arcmotif_dir,eachspe)


        if s0_s6_sub_open_prepare_enrich == 'yes':

            store_prepare_enrich_dir = s0_s6_model_based_motif_enrich_dir + '/store_prepare_enrich_dir'
            if not os.path.exists(store_prepare_enrich_dir):
                os.makedirs(store_prepare_enrich_dir)

            s0_s6_sub_target_spe_str_for_conduct_enrich_list = s0_s6_sub_target_spe_str_for_conduct_enrich_str.split(',')

            for eachspe in s0_s6_sub_target_spe_str_for_conduct_enrich_list:

                store_spe_dir = store_prepare_enrich_dir + '/' + eachspe
                if not os.path.exists(store_spe_dir):
                    os.makedirs(store_spe_dir)


                ipt_acr_motif_modinm_fl = s0_s6_model_based_motif_enrich_dir + '/store_prepare_arcmotif_dir/' + \
                                          '/opt_' + eachspe + '_ACR_motif_modinm.bed'

                ipt_peak_sparse_fl  = input_all_spe_peak_sparse_dir + '/' + 'output_dir_' + eachspe \
                                     + '/opt_peak_sparse_dir/opt_peak_' + eachspe + 'Leaf.sparse'

                s0_subfunctions_motif_model_enrich.subfunction_prepare_regr_enrich_motif(input_required_scripts_dir,
                                                                                         ipt_acr_motif_modinm_fl,
                                                                                         ipt_peak_sparse_fl,
                                                                                         s0_s6_sub_pvalcutoff,
                                                                                         eachspe,store_spe_dir)


        if s0_s6_sub_open_conduct_enrich == 'yes':

            store_conduct_enrich_dir = s0_s6_model_based_motif_enrich_dir + '/store_conduct_enrich_dir'
            if not os.path.exists(store_conduct_enrich_dir):
                os.makedirs(store_conduct_enrich_dir)

            s0_s6_sub_target_spe_str_for_conduct_enrich_list = s0_s6_sub_target_spe_str_for_conduct_enrich_str.split(',')

            for eachspe in s0_s6_sub_target_spe_str_for_conduct_enrich_list:

                opt_spe_dir = store_conduct_enrich_dir + '/' + eachspe
                if not os.path.exists(opt_spe_dir):
                    os.makedirs(opt_spe_dir)

                ipt_meta_add_ACRnum_fl = s0_s6_model_based_motif_enrich_dir + '/store_add_log10tn5_lib_dir/opt_' + eachspe + '.txt'
                ipt_prepare_regr_enrich_dir_dir = s0_s6_model_based_motif_enrich_dir + '/store_prepare_enrich_dir/' + eachspe

                s0_subfunctions_motif_model_enrich.subfunction_conduct_enrichment(input_required_scripts_dir, ipt_meta_add_ACRnum_fl,
                                               ipt_prepare_regr_enrich_dir_dir, opt_spe_dir, s0_s6_sub_target_clusternm,
                                               input_core_num)

    ##updating 121823
    if s0_s7_open_check_cns == 'yes':

        s0_s7_open_check_cns_dir = step00_basic_characters_dir + '/s0_s7_open_check_cns_dir'
        if not os.path.exists(s0_s7_open_check_cns_dir):
            os.makedirs(s0_s7_open_check_cns_dir)

        if s0_s7_sub_open_eachspe_overlapCNS == 'yes':

            store_eachspe_overlapCNS_dir = s0_s7_open_check_cns_dir + '/store_eachspe_overlapCNS_dir'
            if not os.path.exists(store_eachspe_overlapCNS_dir):
                os.makedirs(store_eachspe_overlapCNS_dir)

            ##check the cns
            ##for each species
            s0_s7_target_all_spe_str_list = s0_s7_sub_target_all_spe_str.split(',')

            for eachspe in s0_s7_target_all_spe_str_list:

                store_spe_dir = store_eachspe_overlapCNS_dir + '/' + eachspe
                if not os.path.exists(store_spe_dir):
                    os.makedirs(store_spe_dir)

                ##(ipt_cell_type_spe_fl,ipt_spe_cns_fl, ipt_gene_cate_fl,opt_dir,
                #                       spe_prefix)

                ipt_cell_type_spe_fl = input_all_celltype_acr_dir + '/' + eachspe + '_acr_celltype.txt'
                ipt_spe_cns_fl = input_all_spe_cns_gff_fl_dir + '/ipt_' + eachspe + '_cns_gff.txt'
                ipt_gene_cate_fl = step00_basic_characters_dir + '/s0_s3_open_check_gene_cate_dir/' +eachspe +  '/opt_' + eachspe + '_acr_toGeneCate.txt'
                opt_dir = store_spe_dir
                spe_prefix = eachspe

                s0_subfunctions_cns.subfunction_check_cns_overlap_ACR(ipt_cell_type_spe_fl,ipt_spe_cns_fl, ipt_gene_cate_fl,opt_dir,spe_prefix)

        if s0_s7_sub_open_rice_atlas_CNS == 'yes':

            s0_s7_sub_open_rice_atlas_CNS_dir = s0_s7_open_check_cns_dir + '/s0_s7_sub_open_rice_atlas_CNS_dir'
            if not os.path.exists(s0_s7_sub_open_rice_atlas_CNS_dir):
                os.makedirs(s0_s7_sub_open_rice_atlas_CNS_dir)


            eachspe = 'rice'
            ipt_cell_type_rice_fl = input_all_celltype_acr_dir + '/' + eachspe + '_acr_celltype.txt'
            ipt_atlas_acr_fl = input_total_rice_atlas_acr_fl
            ipt_spe_cns_fl = input_all_spe_cns_gff_fl_dir + '/ipt_' + eachspe + '_cns_gff.txt'
            opt_dir = s0_s7_sub_open_rice_atlas_CNS_dir

            s0_subfunctions_cns.subfunction_intersect_with_atlas(ipt_cell_type_rice_fl, ipt_atlas_acr_fl, ipt_spe_cns_fl, opt_dir)



########
##step01
def step01_species_compare_add_cate (input_rice_to_allspe_blast_dir,
                                     input_all_celltype_acr_dir,
                                     input_all_spe_gene_gff_dir,
                                     input_rice_to_allspe_blast_record_riceID_dir,
                                     input_rice_syntenic_genes_all_os_ACRs_fl,
                                     input_all_syntenic_genes_all_os_ACRs_dir,
                                     input_all_H3K27me3_dir,
                                     input_output_dir,s0_s1_target_spe1,s0_s1_target_spe2_str):

    step01_species_compare_add_cate_dir = input_output_dir + '/step01_species_compare_add_cate_dir'
    if not os.path.exists(step01_species_compare_add_cate_dir):
        os.makedirs(step01_species_compare_add_cate_dir)

    s0_s1_target_spe2_str_list = s0_s1_target_spe2_str.split(',')

    for eachspe2 in s0_s1_target_spe2_str_list:

        #######################
        ##for the rice to maize
        store_rice_to_maize_dir = step01_species_compare_add_cate_dir + '/store_' + s0_s1_target_spe1 + '_to_' + eachspe2 + '_dir'
        if not os.path.exists(store_rice_to_maize_dir):
            os.makedirs(store_rice_to_maize_dir)

        spe1_prefix = s0_s1_target_spe1
        spe2_prefix = eachspe2
        spe1_spe2_syntenic_acr_fl =  input_output_dir + \
                                     '/step00_basic_characters_dir/s0_s1_open_check_syntenic_block_dir' + '/store_' + spe1_prefix + '_to_' + \
                                     spe2_prefix + '_dir/opt_' + spe1_prefix + '_' + spe2_prefix + '_syntenic_region_contain_ACRs.txt'

        input_rice_acr_fl = input_all_celltype_acr_dir + '/' + spe1_prefix + '_acr_celltype.txt'

        if spe1_prefix == 'rice':
            input_rice_to_maize_blast_fl = input_rice_to_allspe_blast_dir + '/' + 'Os_ACRs.vs.' + spe2_prefix + '_blast.bed'
            input_rice_to_maize_blast_record_riceID_fl = input_rice_to_allspe_blast_record_riceID_dir + '/Os_ACRs.vs.' + spe2_prefix + '_blast_record_riceID.bed'

        else:
            input_rice_to_maize_blast_fl = input_rice_to_allspe_blast_dir + '/' + spe1_prefix + '_ACRs.vs.' + 'Os' + '_blast.bed'
            input_rice_to_maize_blast_record_riceID_fl = input_rice_to_allspe_blast_record_riceID_dir + '/' + spe1_prefix + '_ACRs.vs.' + 'Os' + '_blast_record_riceID.bed'


        ##check H3K27me3 if exist
        input_rice_H3K27me3_fl = input_all_H3K27me3_dir + '/' + spe1_prefix + '_H3K27me3_cate.txt'
        input_maize_H3K27me3_fl = input_all_H3K27me3_dir + '/' + spe2_prefix + '_H3K27me3_cate.txt'



        s1_subfunctions.subfunction_iden_correspond(input_rice_to_maize_blast_fl, input_rice_to_maize_blast_record_riceID_fl,
                                    input_rice_H3K27me3_fl, input_maize_H3K27me3_fl,
                                    input_rice_acr_fl, spe1_spe2_syntenic_acr_fl, store_rice_to_maize_dir, spe1_prefix, spe2_prefix)

        input_spe1_gff_fl = input_all_spe_gene_gff_dir + '/' + spe1_prefix + '_gene.gff'
        input_spe2_gff_fl = input_all_spe_gene_gff_dir + '/' + spe2_prefix + '_gene.gff'
        input_summary_fl = store_rice_to_maize_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion.txt'

        input_maize_acr_fl = input_all_celltype_acr_dir + '/' + spe2_prefix + '_acr_celltype.txt'

        input_maize_syntenic_genes_all_os_ACRs_fl = input_all_syntenic_genes_all_os_ACRs_dir + '/' + spe2_prefix + '.syntenic_genes.all_os_ACRs.bed'

        s1_subfunctions.subfunction_check_syntenic_region_gene_pair_direction(input_spe1_gff_fl,input_spe2_gff_fl,
                                                              input_rice_syntenic_genes_all_os_ACRs_fl,input_maize_syntenic_genes_all_os_ACRs_fl,
                                                              input_rice_to_maize_blast_fl,
                                                              input_summary_fl, store_rice_to_maize_dir,input_maize_acr_fl,
                                                              input_rice_H3K27me3_fl,input_maize_H3K27me3_fl,
                                                              spe1_prefix, spe2_prefix)

        ##updating 123123
        ipt_final_summary_fl = store_rice_to_maize_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'
        spe1_H3K27me3_fl = input_rice_H3K27me3_fl
        spe1_acr_fl = input_rice_acr_fl
        s1_subfunctions.subfunction_add_species_specific_one_under_non_syntenic_regions (ipt_final_summary_fl,spe1_H3K27me3_fl,spe1_acr_fl,
                                                                     spe1_prefix,spe2_prefix, store_rice_to_maize_dir)








########
##step02
def step02_summarize_overview(input_output_dir, input_all_syntenic_regions_os_acrs_dir, input_rice_phy_score_fl,
                              input_all_spe_gene_gff_dir, input_all_celltype_acr_dir, input_rice_atlas_acr_coverage_fl,
                              input_total_rice_atlas_acr_fl,
                              input_all_spe_gene_annot_fl_dir, input_all_motif_fimo_dir, input_all_spe_cns_gff_fl_dir,
                              input_marker_gene_list_fl_dir,
                              ipt_cpm_snRNAseq_fl, ipt_maize_cpm_snRNAseq_fl, ipt_cluster_annot_fl,
                              input_all_spe_geneID_symbol_dir, input_all_syntenic_genes_all_os_ACRs_dir,input_all_spe_ortho_dir,input_H3K27me3_fl,
                              input_all_spe_TE_dir,input_DNAmethylation_bw_dir,
                              input_core_num,
                              s0_s1_target_spe1, s0_s1_target_spe2_str,
                              s2_s1_open_basic_summary_overview,
                              s2_s1_open_iden_acr_to_gene_cate,
                              s2_s1_open_check_phy_score,
                              s2_s2_oepn_check_atlasACR_celltypenum,
                              s2_s2_organct_coverage_cutoff,
                              s2_s3_sub_open_only_add_spe1_cpm, s2_s3_open_check_nearby_gene,
                              s2_s3_target_celltype_addCpm_str,
                              s2_s3_sub_open_add_spe1_spe2_cpm, s2_s3_target_spe1_for_twospe_cpm,
                              s2_s3_target_spe2_str_for_twospe_cpm,
                              s2_s3_target_celltype_addCpm_spe2_str,
                              s2_s3_sub_open_add_celltype_specific_str,s2_s3_target_spe1_for_add_ct_str,
                              s2_s3_target_spe2_str_for_add_ct_str,s2_s3_sub_cutoff_avg_fold,
                              s2_s4_open_conduct_enrich_cate_celltype, s2_s4_repeat_times, s2_s4_open_enrich,
                              s2_s5_open_identify_target_motif_binding_plot, s2_s5_target_motif_id_str,
                              s2_s6_open_check_cns_num_per_acr, s2_s6_target_spe2_str, s2_s6_target_spe1,
                              s2_s6_sub_open_check_CNS_overlap_ACR,
                              s2_s6_sub_open_varACRs_in_otherSpes_blastToregionInrice_shown_in_atlas,
                              s2_s7_open_check_snp_num_per_acr, input_rice_3k_snp_fl, s2_s7_target_spe1,
                              s2_s7_target_spe2_str,s2_s8_open_motif_enrich_nonsyn,s2_s8_target_spe1,s2_s8_target_spe2_str,
                              s2_s8_sub_method,s2_s9_open_check_H3K27m3_under_species_specific_ACR,s2_s9_target_spe1,s2_s9_target_spe2_str,
                              s2_s10_open_check_motif_enrich_on_exp_cor_ACR,s2_s10_target_spe1,s2_s10_target_ct_str,s2_s10_repeat_times,
                              s2_s11_open_check_TE_on_ACR,s2_s11_target_spe1,s2_s11_target_spe2_str,
                              s2_s11_sub_open_overlap_enrich_TE,s2_s11_sub_open_motif_enrich_target_TEfam,s2_s11_sub_target_CT_TEfam_str,s2_s11_sub_repeat_times,
                              s2_s11_sub_open_DNAmethylation_target_TEfam
                              ):

    step02_summarize_overview_dir = input_output_dir + '/step02_summarize_overview_dir'
    if not os.path.exists(step02_summarize_overview_dir):
        os.makedirs(step02_summarize_overview_dir)

    if s2_s1_open_basic_summary_overview == 'yes':

        s2_s1_open_basic_summary_overview_dir = step02_summarize_overview_dir + '/s2_s1_open_basic_summary_overview_dir'
        if not os.path.exists(s2_s1_open_basic_summary_overview_dir):
            os.makedirs(s2_s1_open_basic_summary_overview_dir)

        s0_s1_target_spe2_list = s0_s1_target_spe2_str.split(',')

        if s2_s1_open_iden_acr_to_gene_cate == 'yes':

            store_acr_cate_nearbyGene_dir = s2_s1_open_basic_summary_overview_dir + '/store_acr_cate_nearbyGene_dir'
            if not os.path.exists(store_acr_cate_nearbyGene_dir):
                os.makedirs(store_acr_cate_nearbyGene_dir)

            input_rice_acr_fl = input_all_celltype_acr_dir + '/' + s0_s1_target_spe1 + '_acr_celltype.txt'
            input_rice_gff_fl = input_all_spe_gene_gff_dir + '/' + s0_s1_target_spe1 + '_gene.gff'
            s0_subfunctions.subfunction_make_acr_close_to_gene_cate(input_rice_gff_fl, input_rice_acr_fl,
                                                                    store_acr_cate_nearbyGene_dir, s0_s1_target_spe1)

            for eachspe2 in s0_s1_target_spe2_list:

                input_spe2_acr_fl = input_all_celltype_acr_dir + '/' + eachspe2 + '_acr_celltype.txt'
                input_spe2_gff_fl =  input_all_spe_gene_gff_dir + '/' + eachspe2 + '_gene.gff'

                s0_subfunctions.subfunction_make_acr_close_to_gene_cate(input_spe2_gff_fl, input_spe2_acr_fl,
                                                                        store_acr_cate_nearbyGene_dir, eachspe2)


        for eachspe2 in s0_s1_target_spe2_list:

            #############################
            ##for the rice to maize CTver
            store_rice_to_maize_CTver_dir = s2_s1_open_basic_summary_overview_dir + '/store_' + s0_s1_target_spe1 + '_to_' + eachspe2 + '_CTver_dir'
            if not os.path.exists(store_rice_to_maize_CTver_dir):
                os.makedirs(store_rice_to_maize_CTver_dir)

            spe1_prefix = s0_s1_target_spe1
            spe2_prefix = eachspe2

            ipt_final_summary_blast_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                         spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'
            ipt_sep2_syntenic_region_fl = input_all_syntenic_regions_os_acrs_dir + '/' + spe2_prefix + '.syntenic_regions.os_acrs.bed'
            opt_dir = store_rice_to_maize_CTver_dir

            ipt_spe1_acr_cate_fl = s2_s1_open_basic_summary_overview_dir + '/store_acr_cate_nearbyGene_dir/opt_' + spe1_prefix + '_acr_toGeneCate.txt'
            ipt_spe2_acr_cate_fl = s2_s1_open_basic_summary_overview_dir + '/store_acr_cate_nearbyGene_dir/opt_' + spe2_prefix + '_acr_toGeneCate.txt'

            input_spe1_acr_fl = input_all_celltype_acr_dir + '/' + spe1_prefix + '_acr_celltype.txt'

            s2_decide_celltype_cate = 'CTVer'
            s2_subfunctions.subfunction_summarize_overview(ipt_final_summary_blast_fl, ipt_sep2_syntenic_region_fl, input_rice_phy_score_fl,opt_dir,
                                                           ipt_spe1_acr_cate_fl, ipt_spe2_acr_cate_fl,input_spe1_acr_fl,
                                                           spe1_prefix,spe2_prefix,s2_decide_celltype_cate,s2_s1_open_check_phy_score)

            ipt_sharednotsharedCate_celltype_fl = s2_s1_open_basic_summary_overview_dir + '/' + 'store_' + spe1_prefix + '_to_' + spe2_prefix + '_CTver_dir/opt1_syntenic_region_' + \
                                                  spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_celltypes.txt'
            s2_subfunctions.subfunction_check_not_shared_statistic_different_enrich (ipt_sharednotsharedCate_celltype_fl,opt_dir,spe1_prefix,spe2_prefix)


    if s2_s2_oepn_check_atlasACR_celltypenum == 'yes':

        s2_s2_oepn_check_atlasACR_celltypenum_dir = step02_summarize_overview_dir + '/s2_s2_oepn_check_atlasACR_celltypenum_dir'
        if not os.path.exists(s2_s2_oepn_check_atlasACR_celltypenum_dir):
            os.makedirs(s2_s2_oepn_check_atlasACR_celltypenum_dir)

        s0_s1_target_spe2_str_list = s0_s1_target_spe2_str.split(',')

        for eachspe2 in s0_s1_target_spe2_str_list:


            spe1_prefix = s0_s1_target_spe1
            spe2_prefix = eachspe2

            store_opt_dir = s2_s2_oepn_check_atlasACR_celltypenum_dir + '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir'
            if not os.path.exists(store_opt_dir):
                os.makedirs(store_opt_dir)


            ipt_rice_leaf_syntenic_fl = step02_summarize_overview_dir +  '/s2_s1_open_basic_summary_overview_dir' + \
                                        '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_CTver_dir' + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + '_SharedNotSharedCate_ACRs_sorted.txt'

            #ipt_rice_leaf_syntenic_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + \
            #                            spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
            #                            spe1_prefix + '_' + spe2_prefix+ '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

            ipt_rice_leaf_syntenic_addCTcate_fl = step02_summarize_overview_dir + '/s2_s1_open_basic_summary_overview_dir' + \
                                        '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_CTver_dir' + '/opt1_syntenic_region_' + spe1_prefix + '_' + spe2_prefix + \
                                                  '_SharedNotSharedCate_ACRs_addCTcate_sorted.txt'

            s2_subfunctions.subfunction_check_acc_in_rice_atlas_peaks(input_rice_atlas_acr_coverage_fl, input_total_rice_atlas_acr_fl,
                                                      ipt_rice_leaf_syntenic_fl,ipt_rice_leaf_syntenic_addCTcate_fl,
                                                      store_opt_dir,
                                                      spe1_prefix, spe2_prefix,
                                                      s2_s2_organct_coverage_cutoff)


    if s2_s3_open_check_nearby_gene == 'yes':

        s2_s3_open_check_nearby_gene_dir = step02_summarize_overview_dir + '/s2_s3_open_check_nearby_gene_dir'
        if not os.path.exists(s2_s3_open_check_nearby_gene_dir):
            os.makedirs(s2_s3_open_check_nearby_gene_dir)

        ##we will focuse on the rice gene
        ##for the shared-Acc
        ##what is the rice gene cloest to the broad shared-Acc
        ##what is the rice gene cloese to the

        ##how many of the genes are identified in all five species
        spe1_prefix = s0_s1_target_spe1


        s0_s1_target_spe2_str_list = s0_s1_target_spe2_str.split(',')

        if s2_s3_sub_open_only_add_spe1_cpm == 'yes':

            for eachspe2 in s0_s1_target_spe2_str_list:

                spe2_prefix = eachspe2

                opt_store_spe_dir = s2_s3_open_check_nearby_gene_dir + '/' + eachspe2
                if not os.path.exists(opt_store_spe_dir):
                    os.makedirs(opt_store_spe_dir)

                ipt_final_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                             spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype_addSpeSpec.txt'

                ipt_spe1_gff_fl = input_all_spe_gene_gff_dir + '/' + spe1_prefix + '_gene.gff'
                ipt_spe1_acr_sumit_gene_closest_fl = step02_summarize_overview_dir + '/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/' + \
                                                    'temp_' + spe1_prefix + '_acr_sumit_gene_closest.txt'

                ipt_spe1_gene_annot_fl = input_all_spe_gene_annot_fl_dir + '/' + spe1_prefix + '_gene_annot.txt'

                s2_subfunctions.subfunction_check_nearby_genes (ipt_final_summary_fl,ipt_spe1_gff_fl,ipt_spe1_acr_sumit_gene_closest_fl,
                                                                 ipt_spe1_gene_annot_fl,opt_store_spe_dir,spe1_prefix,spe2_prefix)

                ipt_marker_gene_fl = input_marker_gene_list_fl_dir + '/' + spe1_prefix + '_markers.txt'
                ipt_final_add_gene_fl = opt_store_spe_dir + '/opt_' + spe1_prefix + '_'  + spe2_prefix + '_final_blast_summary_add_gene.txt'
                s2_subfunctions.subfunction_add_marker_gene_information(ipt_marker_gene_fl,ipt_final_add_gene_fl,opt_store_spe_dir,spe1_prefix,spe2_prefix)

                ipt_final_add_gene_fl = opt_store_spe_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene.txt'
                s2_subfunctions.subfunction_add_gene_exp_data_per_celltype(ipt_cpm_snRNAseq_fl, ipt_cluster_annot_fl,
                                                                           ipt_final_add_gene_fl, s2_s3_target_celltype_addCpm_str,
                                                                           opt_store_spe_dir,
                                                                            spe1_prefix, spe2_prefix)

                ipt_summary_celltypeCpm_fl = opt_store_spe_dir + '/' + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_add_celltypeCpm.txt'
                ipt_spe1_acr_cate_fl = step02_summarize_overview_dir  + '/s2_s1_open_basic_summary_overview_dir' + '/store_acr_cate_nearbyGene_dir/opt_' + spe1_prefix + '_acr_toGeneCate.txt'
                s2_subfunctions.subfunction_add_distance_to_celltypeCpm(ipt_summary_celltypeCpm_fl, opt_store_spe_dir, ipt_spe1_acr_cate_fl,
                                                        spe1_prefix, spe2_prefix)


        if s2_s3_sub_open_add_spe1_spe2_cpm == 'yes':

            spe1_prefix = s2_s3_target_spe1_for_twospe_cpm

            for eachspe2 in s2_s3_target_spe2_str_for_twospe_cpm.split(','):

                spe2_prefix = eachspe2

                opt_store_spe_dir = s2_s3_open_check_nearby_gene_dir + '/' + eachspe2
                if not os.path.exists(opt_store_spe_dir):
                    os.makedirs(opt_store_spe_dir)


                ##updating 010924 we will incorporate the non-syntenic regions
                #ipt_final_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                #                             spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'


                ipt_final_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                             spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype_addSpeSpec.txt'



                ipt_spe1_gff_fl = input_all_spe_gene_gff_dir + '/' + spe1_prefix + '_gene.gff'
                ipt_spe1_acr_sumit_gene_closest_fl = step02_summarize_overview_dir + '/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/' + \
                                                     'temp_' + spe1_prefix + '_acr_sumit_gene_closest.txt'
                ipt_spe1_gene_annot_fl = input_all_spe_gene_annot_fl_dir + '/' + spe1_prefix + '_gene_annot.txt'

                ipt_spe2_gff_fl = input_all_spe_gene_gff_dir + '/' + spe2_prefix + '_gene.gff'
                ipt_spe1_geneID_symbol_fl = input_all_spe_geneID_symbol_dir + '/' + 'opt_' + spe1_prefix + '_symbol.txt'

                store_region_sumit_dir = opt_store_spe_dir + '/' + 'store_region_sumit_dir'
                if not os.path.exists(store_region_sumit_dir):
                    os.makedirs(store_region_sumit_dir)

                s2_subfunctions.subfunction_make_region_close_to_gene_cate(ipt_spe2_gff_fl, ipt_final_summary_fl, store_region_sumit_dir, spe2_prefix)

                ipt_spe2_region_gene_closest_fl =  store_region_sumit_dir + '/temp_' + spe2_prefix + '_region_sumit_gene_closest.txt'
                s2_subfunctions.subfunction_check_nearby_genes_two_spe (ipt_final_summary_fl,ipt_spe1_gff_fl,ipt_spe1_acr_sumit_gene_closest_fl,ipt_spe1_gene_annot_fl,
                                            ipt_spe2_gff_fl,ipt_spe2_region_gene_closest_fl,ipt_spe1_geneID_symbol_fl,
                                            opt_store_spe_dir,spe1_prefix,spe2_prefix)

                ipt_final_summary_with_geneID_two_spe_fl = opt_store_spe_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene_two_Spe.txt'
                s2_subfunctions.subfunction_add_gene_exp_data_per_celltype_two_spe(ipt_cpm_snRNAseq_fl, ipt_cluster_annot_fl,
                                                                                   ipt_maize_cpm_snRNAseq_fl,
                                                                   ipt_final_summary_with_geneID_two_spe_fl,
                                                                   s2_s3_target_celltype_addCpm_str,
                                                                   s2_s3_target_celltype_addCpm_spe2_str, opt_store_spe_dir,
                                                                   spe1_prefix, spe2_prefix)

                #ipt_final_summary_add_gene_ID_fl = opt_store_spe_dir + '/opt_' + spe1_prefix + '_'  + spe2_prefix + '_final_blast_summary_add_gene.txt'
                #ipt_spe2_syntenic_region_gene_fl = input_all_syntenic_genes_all_os_ACRs_dir + '/' + spe2_prefix + '.syntenic_genes.all_os_ACRs.bed'

                ##we will not analyze this as the ortho is not the right one
                #ipt_ortho_fl = input_all_spe_ortho_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '.txt'

                #s2_subfunctions.subfunction_add_gene_exp_data_per_celltype_two_spe_species_specificVer (ipt_final_summary_add_gene_ID_fl,ipt_cluster_annot_fl,ipt_cpm_snRNAseq_spe1_fl,
                #                                                            ipt_cpm_snRNAseq_spe2_fl,ipt_ortho_fl,
                #                                                            ipt_spe2_syntenic_region_gene_fl,opt_store_spe_dir,
                #                                                            s2_s3_target_celltype_addCpm_str,s2_s3_target_celltype_addCpm_spe2_str,
                #                                                            spe1_prefix,spe2_prefix)


        ##udpating 010924
        ##we will add the cell type specific str to the gene file
        if s2_s3_sub_open_add_celltype_specific_str == 'yes':

            spe1_prefix = s2_s3_target_spe1_for_add_ct_str

            for eachspe2 in s2_s3_target_spe2_str_for_add_ct_str.split(','):

                spe2_prefix = eachspe2

                opt_store_spe_dir = s2_s3_open_check_nearby_gene_dir + '/' + eachspe2
                if not os.path.exists(opt_store_spe_dir):
                    os.makedirs(opt_store_spe_dir)

                ipt_spe1_spe2_orth_fl = input_all_spe_ortho_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '.txt'
                opt_dir = opt_store_spe_dir
                s2_s3_target_celltype_addCpm_spe1_str = s2_s3_target_celltype_addCpm_str

                ##check the ortholog
                if spe2_prefix == 'maize':

                    ipt_cpm_snRNAseq_spe1_fl = ipt_cpm_snRNAseq_fl
                    ipt_cpm_snRNAseq_spe2_fl = ipt_maize_cpm_snRNAseq_fl

                else:
                    ipt_cpm_snRNAseq_spe1_fl = ipt_maize_cpm_snRNAseq_fl
                    ipt_cpm_snRNAseq_spe2_fl = ipt_cpm_snRNAseq_fl




                s2_subfunctions_exp.subfunction_check_ortho_non_syn_synSS_class_gene_nearby_ACRs (ipt_spe1_spe2_orth_fl,ipt_cluster_annot_fl,
                                                                                                  ipt_cpm_snRNAseq_spe1_fl, ipt_cpm_snRNAseq_spe2_fl,opt_dir,
                                                                  s2_s3_target_celltype_addCpm_spe1_str,s2_s3_target_celltype_addCpm_spe2_str,s2_s3_sub_cutoff_avg_fold,
                                                                  spe1_prefix,spe2_prefix)

                ##add the cloest gene information
                ipt_final_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                       spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype_addSpeSpec.txt'

                ipt_spe1_gff_fl = input_all_spe_gene_gff_dir + '/' + spe1_prefix + '_gene.gff'
                ipt_spe1_acr_sumit_gene_closest_fl = step02_summarize_overview_dir + '/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/' + \
                                                     'temp_' + spe1_prefix + '_acr_sumit_gene_closest.txt'

                ipt_spe1_gene_annot_fl = input_all_spe_gene_annot_fl_dir + '/' + spe1_prefix + '_gene_annot.txt'

                s2_subfunctions.subfunction_check_nearby_genes(ipt_final_summary_fl, ipt_spe1_gff_fl,
                                                               ipt_spe1_acr_sumit_gene_closest_fl,
                                                               ipt_spe1_gene_annot_fl, opt_store_spe_dir, spe1_prefix,
                                                               spe2_prefix)

                ##add the ct specific
                ipt_final_summary_fl = opt_dir + '/opt_' + spe1_prefix + '_' + spe2_prefix + '_final_blast_summary_add_gene.txt'
                ipt_ortho_gene_celltype_specific_fl = opt_dir + '/opt_ortho_gene_celltype_specific_exp' + spe1_prefix + '_' + spe2_prefix + '.txt'

                s2_subfunctions_exp.subfunction_add_orth_celltype_to_final_summary_fl(ipt_final_summary_fl,
                                                                                      ipt_ortho_gene_celltype_specific_fl, spe1_prefix,spe2_prefix, opt_dir)


    if s2_s4_open_conduct_enrich_cate_celltype == 'yes':

        s2_s4_open_conduct_enrich_cate_celltype_dir = step02_summarize_overview_dir + '/s2_s4_open_conduct_enrich_cate_celltype_dir'
        if not os.path.exists(s2_s4_open_conduct_enrich_cate_celltype_dir):
            os.makedirs(s2_s4_open_conduct_enrich_cate_celltype_dir)

        s0_s1_target_spe2_str_list = s0_s1_target_spe2_str.split(',')

        for eachspe2 in s0_s1_target_spe2_str_list:

            spe1_prefix = s0_s1_target_spe1
            spe2_prefix = eachspe2

            store_opt_dir = s2_s4_open_conduct_enrich_cate_celltype_dir + '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir'
            if not os.path.exists(store_opt_dir):
                os.makedirs(store_opt_dir)

            store_opt_dir_for_sharedA = s2_s4_open_conduct_enrich_cate_celltype_dir + '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_for_sharedA_dir'
            if not os.path.exists(store_opt_dir_for_sharedA):
                os.makedirs(store_opt_dir_for_sharedA)

            ipt_total_spe_acr_fl = input_all_celltype_acr_dir + '/' + spe1_prefix + '_acr_celltype.txt'

            ipt_motif_fl = input_all_motif_fimo_dir + '/' + spe1_prefix + '_motif_fimo.txt'

            ipt_final_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                         spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

            if s2_s4_open_enrich == 'yes':

               #s2_subfunctions.subfunction_check_motif_enrichment_per_celltype (ipt_motif_fl,ipt_final_summary_fl,ipt_total_spe_acr_fl,store_opt_dir,
               #                                          s2_s4_repeat_times,input_core_num)

                s2_subfunctions.subfunction_check_motif_enrichment_per_celltype_sharedA(ipt_motif_fl, ipt_final_summary_fl,
                                                                        ipt_total_spe_acr_fl, store_opt_dir_for_sharedA,
                                                                        s2_s4_repeat_times, input_core_num)



            ipt_spe_dir = store_opt_dir
            opt_dir = s2_s4_open_conduct_enrich_cate_celltype_dir
            s2_subfunctions.subfunction_collect_all_enrich_results(ipt_spe_dir, opt_dir, spe1_prefix, spe2_prefix)



            ipt_spe_dir = store_opt_dir_for_sharedA
            opt_dir = s2_s4_open_conduct_enrich_cate_celltype_dir
            s2_subfunctions.subfunction_collect_all_enrich_results(ipt_spe_dir, opt_dir, spe1_prefix, spe2_prefix)


    if s2_s5_open_identify_target_motif_binding_plot == 'yes':

        s2_s5_open_identify_target_motif_binding_plot = step02_summarize_overview_dir + '/s2_s5_open_identify_target_motif_binding_plot'
        if not os.path.exists(s2_s5_open_identify_target_motif_binding_plot):
            os.makedirs(s2_s5_open_identify_target_motif_binding_plot)

        store_final_summary_all_spe_fl_dir = s2_s5_open_identify_target_motif_binding_plot + '/store_final_summary_all_spe_fl_dir'
        if not os.path.exists(store_final_summary_all_spe_fl_dir):
            os.makedirs(store_final_summary_all_spe_fl_dir)


        s0_s1_target_spe2_str_list = s0_s1_target_spe2_str.split(',')

        for eachspe2 in s0_s1_target_spe2_str_list:

            spe1_prefix = s0_s1_target_spe1
            spe2_prefix = eachspe2

            ipt_final_summary_blast_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                         spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

            cmd = 'cp ' + ipt_final_summary_blast_fl + ' ' + store_final_summary_all_spe_fl_dir + '/' + spe2_prefix + '_summary.txt'
            print(cmd)
            subprocess.call(cmd,shell=True)

        ipt_motif_spe1_fimo_fl = input_all_motif_fimo_dir + '/' + s0_s1_target_spe1 + '_motif_fimo.txt'

        s2_subfunctions.subfunction_identify_target_motif_binding_plotting (store_final_summary_all_spe_fl_dir,
                                                                           ipt_motif_spe1_fimo_fl, s2_s5_target_motif_id_str,
                                                                            s2_s5_open_identify_target_motif_binding_plot,
                                                                            s0_s1_target_spe1)


    if s2_s6_open_check_cns_num_per_acr == 'yes':

        s2_s6_open_check_cns_num_per_acr = step02_summarize_overview_dir + '/s2_s6_open_check_cns_num_per_acr_dir'
        if not os.path.exists(s2_s6_open_check_cns_num_per_acr):
            os.makedirs(s2_s6_open_check_cns_num_per_acr)

        s2_s5_target_spe2_str_list = s2_s6_target_spe2_str.split(',')

        for eachspe2 in s2_s5_target_spe2_str_list:

            spe1_prefix = s2_s6_target_spe1
            spe2_prefix = eachspe2

            store_opt_dir = s2_s6_open_check_cns_num_per_acr + '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir'
            if not os.path.exists(store_opt_dir):
                os.makedirs(store_opt_dir)

            ipt_spe1_cns_gff_fl = input_all_spe_cns_gff_fl_dir + '/ipt_' + spe1_prefix + '_cns_gff.txt'
            ipt_final_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                         spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype_addSpeSpec.txt'

            ##updating 121923
            if s2_s6_sub_open_check_CNS_overlap_ACR == 'yes':


                ##updating 020224
                ##we will not use the species specific cases
                ipt_not_speSpec_final_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                         spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

                s2_subfunctions_cns.subfunction_check_CNS_overlap_density(ipt_spe1_cns_gff_fl, ipt_not_speSpec_final_summary_fl, store_opt_dir)

                ipt_final_cate_acr_cns_num_fl = store_opt_dir + '/opt_final_cate_acr_cns_num.txt'
                s2_subfunctions_cns.subfunction_check_enrichment_for_celltype (ipt_final_cate_acr_cns_num_fl,store_opt_dir)




            if s2_s6_sub_open_varACRs_in_otherSpes_blastToregionInrice_shown_in_atlas == 'yes':

                ipt_rice_atlas_fl = input_total_rice_atlas_acr_fl
                opt_dir = store_opt_dir
                ipt_spe2_cns_gff_fl = input_all_spe_cns_gff_fl_dir + '/ipt_' + spe2_prefix + '_cns_gff.txt'

                ipt_not_speSpec_final_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                                   spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

                s2_subfunctions_cns.subfunction_check_varACRotherSpe_inRiceRegion_in_Atlas (ipt_not_speSpec_final_summary_fl,ipt_rice_atlas_fl,ipt_spe2_cns_gff_fl,input_rice_atlas_acr_coverage_fl,
                                                            opt_dir,spe1_prefix,spe2_prefix,s2_s2_organct_coverage_cutoff)

                ipt_summary_fl = store_opt_dir + '/opt_' + spe1_prefix + '_to_' + spe2_prefix + '_summary_add_' + spe2_prefix + '_atlas_acr.txt'

                s2_subfunctions_cns.subfunction_add_check_num_celltypeACR(ipt_summary_fl, ipt_rice_atlas_fl,
                                                      input_rice_atlas_acr_coverage_fl,
                                                      s2_s2_organct_coverage_cutoff, spe1_prefix, spe2_prefix, opt_dir)



    ##updating 122723
    if s2_s7_open_check_snp_num_per_acr == 'yes':

        s2_s7_open_check_snp_num_per_acr_dir = step02_summarize_overview_dir + '/s2_s7_open_check_snp_num_per_acr_dir'
        if not os.path.exists(s2_s7_open_check_snp_num_per_acr_dir):
            os.makedirs(s2_s7_open_check_snp_num_per_acr_dir)

        s2_s7_target_spe2_str_list = s2_s7_target_spe2_str.split(',')

        for eachspe2 in s2_s7_target_spe2_str_list:

            spe1_prefix = s2_s7_target_spe1
            spe2_prefix = eachspe2

            store_opt_dir = s2_s7_open_check_snp_num_per_acr_dir + '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir'
            if not os.path.exists(store_opt_dir):
                os.makedirs(store_opt_dir)

            ipt_snp_fl = input_rice_3k_snp_fl
            ipt_final_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                         spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'
            opt_dir = store_opt_dir

            s2_subfunctions_snp.subfunction_add_snp_num_to_acr (ipt_snp_fl,ipt_final_summary_fl,spe1_prefix,spe2_prefix,opt_dir)


    if s2_s8_open_motif_enrich_nonsyn == 'yes':

        s2_s8_open_motif_enrich_nonsyn_dir = step02_summarize_overview_dir + '/s2_s8_open_motif_enrich_nonsyn_dir'
        if not os.path.exists(s2_s8_open_motif_enrich_nonsyn_dir):
            os.makedirs(s2_s8_open_motif_enrich_nonsyn_dir)

        s2_s8_target_spe2_str_list = s2_s8_target_spe2_str.split(',')

        for eachspe2 in s2_s8_target_spe2_str_list:

            spe1_prefix = s2_s8_target_spe1
            spe2_prefix = eachspe2

            store_opt_dir = s2_s8_open_motif_enrich_nonsyn_dir + '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir'
            if not os.path.exists(store_opt_dir):
                os.makedirs(store_opt_dir)

            ipt_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir' + \
                             '/store_' + spe1_prefix + '_to_' + eachspe2 + '_dir/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype_addSpeSpec.txt'

            ipt_fimo_motif_prediction_fl = input_all_motif_fimo_dir + '/' + spe1_prefix + '_motif_fimo.txt'


            if s2_s8_sub_method == 'fisher':

                store_fisher_results_dir = store_opt_dir + '/store_fisher_results_dir'
                if not os.path.exists(store_fisher_results_dir):
                    os.makedirs(store_fisher_results_dir)

                s2_subfunctions_nonsyn_syn.subfunction_fisher_motif_enrichment_analysis_for_syn_nonsyn (ipt_fimo_motif_prediction_fl, ipt_summary_fl,
                                                                    spe1_prefix, spe2_prefix, store_fisher_results_dir)

            if s2_s8_sub_method == 'binomial':
                ipt_total_spe_acr_fl = input_all_celltype_acr_dir + '/' + spe1_prefix + '_acr_celltype.txt'

                store_binomial_results_dir = store_opt_dir + '/store_binomial_results_dir'
                if not os.path.exists(store_binomial_results_dir):
                    os.makedirs(store_binomial_results_dir)

                s2_subfunctions_nonsyn_syn.subfunction_check_motif_enrichment_per_celltype_binomial(ipt_fimo_motif_prediction_fl,ipt_summary_fl,ipt_total_spe_acr_fl,store_binomial_results_dir,s2_s4_repeat_times,input_core_num)

                s2_subfunctions_nonsyn_syn.subfunction_collect_enrich_binomial_results (store_binomial_results_dir,store_opt_dir,spe1_prefix,spe2_prefix)

    ##udpating 010624
    if s2_s9_open_check_H3K27m3_under_species_specific_ACR == 'yes':

        s2_s9_open_check_H3K27m3_under_species_specific_ACR_dir = step02_summarize_overview_dir + '/s2_s9_open_check_H3K27m3_under_species_specific_ACR_dir'
        if not os.path.exists(s2_s9_open_check_H3K27m3_under_species_specific_ACR_dir):
            os.makedirs(s2_s9_open_check_H3K27m3_under_species_specific_ACR_dir)

        s2_s9_target_spe2_list = s2_s9_target_spe2_str.split(',')

        store_spe_ipt_summary_fl_dic = {}
        for eachspe2 in s2_s9_target_spe2_list:

            spe1_prefix = s2_s9_target_spe1
            spe2_prefix = eachspe2

            store_opt_dir = s2_s9_open_check_H3K27m3_under_species_specific_ACR_dir + '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir'
            if not os.path.exists(store_opt_dir):
                os.makedirs(store_opt_dir)

            ipt_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir' + \
                             '/store_' + spe1_prefix + '_to_' + eachspe2 + '_dir/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype_addSpeSpec.txt'

            store_spe_ipt_summary_fl_dic[eachspe2] = ipt_summary_fl

        s2_subfunctions_H3K27me3.subfunction_plot_H3K27me3_on_species_specific_ACRs(store_spe_ipt_summary_fl_dic, input_H3K27me3_fl, s2_s9_open_check_H3K27m3_under_species_specific_ACR_dir,
                                                           input_core_num)

    ##udpating 011024
    if s2_s10_open_check_motif_enrich_on_exp_cor_ACR == 'yes':

        s2_s10_open_check_motif_enrich_on_exp_cor_ACR_dir = step02_summarize_overview_dir + '/s2_s10_open_check_motif_enrich_on_exp_cor_ACR_dir'
        if not os.path.exists(s2_s10_open_check_motif_enrich_on_exp_cor_ACR_dir):
            os.makedirs(s2_s10_open_check_motif_enrich_on_exp_cor_ACR_dir)

        s2_s9_target_spe2_list = s2_s9_target_spe2_str.split(',')


        for eachspe2 in s2_s9_target_spe2_list:

            spe1_prefix = s2_s10_target_spe1
            spe2_prefix = eachspe2

            store_opt_dir = s2_s10_open_check_motif_enrich_on_exp_cor_ACR_dir + '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir'
            if not os.path.exists(store_opt_dir):
                os.makedirs(store_opt_dir)

            ipt_final_summary_fl = step02_summarize_overview_dir + '/s2_s3_open_check_nearby_gene_dir/' + eachspe2 + \
                                   '/opt_' + spe1_prefix + '_' + eachspe2 + '_final_blast_summary_add_gene_add_synOrNsyn_orth_exp_CTstr_collapes.txt'

            store_target_acr_list_fl_dir = store_opt_dir + '/' + 'store_target_acr_list_fl_dir'
            if not os.path.exists(store_target_acr_list_fl_dir):
                os.makedirs(store_target_acr_list_fl_dir)


            s2_subfunctions_exp.subfunction_prepare_target_acr_loc(ipt_final_summary_fl, store_target_acr_list_fl_dir, s2_s10_target_ct_str,spe1_prefix)

            ipt_target_celltype_acr_loc_fl_list = glob.glob(store_target_acr_list_fl_dir + '/opt_*')

            ipt_total_spe_acr_fl = input_all_celltype_acr_dir + '/' + spe1_prefix + '_acr_celltype.txt'
            ipt_motif_fl = input_all_motif_fimo_dir + '/' + spe1_prefix + '_motif_fimo.txt'

            s2_subfunctions_exp.subfunction_check_motif_enrichment(ipt_target_celltype_acr_loc_fl_list, ipt_total_spe_acr_fl,
                                               ipt_motif_fl, s2_s10_repeat_times, store_opt_dir)



    ##updating 011124
    if s2_s11_open_check_TE_on_ACR == 'yes':

        s2_s11_open_check_TE_on_ACR_dir = step02_summarize_overview_dir + '/s2_s11_open_check_TE_on_ACR_dir'
        if not os.path.exists(s2_s11_open_check_TE_on_ACR_dir):
            os.makedirs(s2_s11_open_check_TE_on_ACR_dir)

        s2_s11_target_spe2_list = s2_s11_target_spe2_str.split(',')

        for eachspe2 in s2_s11_target_spe2_list:
            spe1_prefix = s2_s11_target_spe1
            spe2_prefix = eachspe2

            store_opt_dir = s2_s11_open_check_TE_on_ACR_dir + '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir'
            if not os.path.exists(store_opt_dir):
                os.makedirs(store_opt_dir)

            ipt_final_summary_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + \
                                         spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype_addSpeSpec.txt'

            ipt_spe1_TE_fl = input_all_spe_TE_dir + '/' + spe1_prefix + '_TE.gff3'

            if s2_s11_sub_open_overlap_enrich_TE == 'yes':

                s2_subfunctions_nonsyn_syn.subfunction_check_TE_all_fam_ACR_composition_syn_nonsyn(ipt_spe1_TE_fl, ipt_final_summary_fl, store_opt_dir)

                ipt_count_prop_fl = store_opt_dir + '/opt_final_TE_prop_in_eachcate.txt'
                s2_subfunctions_nonsyn_syn.subfunction_conduct_enrichment_from_prop_fl_fisher (ipt_count_prop_fl,store_opt_dir)

                ##for the cell type enrichment
                s2_subfunctions_nonsyn_syn.subfunction_check_TE_all_fam_ACR_composition_syn_nonsyn_per_celltype(ipt_spe1_TE_fl, ipt_final_summary_fl,
                                                                                                                store_opt_dir)
                ipt_count_prop_fl = store_opt_dir + '/opt_final_TE_prop_in_eachcate_celltypeVer.txt'
                s2_subfunctions_nonsyn_syn.subfunction_conduct_enrichment_from_prop_fl_fisher_per_celltype(ipt_count_prop_fl, store_opt_dir)


            if s2_s11_sub_open_motif_enrich_target_TEfam == 'yes':

                store_target_TE_family_enrich_motif_dir = store_opt_dir + '/store_target_TE_family_enrich_motif_dir'
                if not os.path.exists(store_target_TE_family_enrich_motif_dir):
                    os.makedirs(store_target_TE_family_enrich_motif_dir)

                ipt_target_CT_TEfam_str = s2_s11_sub_target_CT_TEfam_str
                ipt_total_spe_acr_fl = input_all_celltype_acr_dir + '/' + spe1_prefix + '_acr_celltype.txt'
                ipt_motif_fl = input_all_motif_fimo_dir + '/' + spe1_prefix + '_motif_fimo.txt'

                s2_subfunctions_nonsyn_syn.subfunction_conduct_motif_enrichment_target_CT_TEfam(ipt_target_CT_TEfam_str,
                                                                                                store_opt_dir, ipt_total_spe_acr_fl,
                                                                     ipt_motif_fl, store_target_TE_family_enrich_motif_dir,
                                                                     s2_s11_sub_repeat_times)

            if s2_s11_sub_open_DNAmethylation_target_TEfam == 'yes':

                store_target_TE_family_DNAmethy_dir = store_opt_dir + '/store_target_TE_family_DNAmethy_dir'
                if not os.path.exists(store_target_TE_family_DNAmethy_dir):
                    os.makedirs(store_target_TE_family_DNAmethy_dir)

                input_all_three_methy_context_dir = input_DNAmethylation_bw_dir
                ipt_target_CT_TEfam_str = s2_s11_sub_target_CT_TEfam_str
                opt_dir_store_intersect_dir = store_opt_dir
                ipt_total_spe_acr_fl = input_all_celltype_acr_dir + '/' + spe1_prefix + '_acr_celltype.txt'
                opt_dir = store_target_TE_family_DNAmethy_dir
                s2_subfunctions_nonsyn_syn.subfunction_conduct_DNA_methylation_target_CT_TEfam(input_all_three_methy_context_dir,
                                                                    ipt_target_CT_TEfam_str,
                                                                    opt_dir_store_intersect_dir,
                                                                    ipt_spe1_TE_fl, ipt_total_spe_acr_fl,
                                                                    input_core_num, opt_dir)



########
##step03
def step03_summarize_overview_H3K27me3 (input_output_dir,input_all_celltype_acr_dir,input_all_spe_gene_gff_dir,
                                        input_target_organ_H3K27_peak_all_dir,input_rice_all_organ_H3K27_fl_dir,
                                        input_all_spe_PRE_motif_fimo_dir,input_all_spe_cns_gff_fl_dir,
                                        input_total_rice_atlas_acr_fl,
                                        s3_s1_open_iden_acr_to_gene_cate,
                                        s3_s2_open_summarize_overview_cross_species,
                                        s3_s2_extend_H3K27me3_flank_bp,s3_s2_target_spe1,s3_s2_target_spe2_str,
                                        s3_s3_open_check_conserved_H3K27me3_broad_ACR_in_other_tissue,
                                        s3_s3_target_organ_str,
                                        s3_s4_open_check_motif_freq_in_diff_cate_acr,
                                        s3_s5_open_check_cns_overlap_ACR,s3_s6_open_check_overlap_atlas_acr,
                                        s3_s7_open_check_target_motif_in_ACR,s3_s7_target_motifID):

    step03_summarize_overview_H3K27me3_dir = input_output_dir + '/step03_summarize_overview_H3K27me3_dir'
    if not os.path.exists(step03_summarize_overview_H3K27me3_dir):
        os.makedirs(step03_summarize_overview_H3K27me3_dir)


    if s3_s1_open_iden_acr_to_gene_cate == 'yes':

        s3_1_store_acr_cate_nearbyGene_dir = step03_summarize_overview_H3K27me3_dir + '/s3_1_store_acr_cate_nearbyGene_dir'
        if not os.path.exists(s3_1_store_acr_cate_nearbyGene_dir):
            os.makedirs(s3_1_store_acr_cate_nearbyGene_dir)

        input_rice_acr_fl = input_all_celltype_acr_dir + '/' + 'rice' + '_acr_celltype.txt'
        input_maize_acr_fl = input_all_celltype_acr_dir + '/' + 'maize' + '_acr_celltype.txt'
        input_sorghum_acr_fl = input_all_celltype_acr_dir + '/' + 'sorghum' + '_acr_celltype.txt'

        input_rice_gff_fl = input_all_spe_gene_gff_dir + '/' + 'rice' + '_gene.gff'
        input_maize_gff_fl = input_all_spe_gene_gff_dir + '/' + 'maize' + '_gene.gff'
        input_sorghum_gff_fl = input_all_spe_gene_gff_dir + '/' + 'sorghum' + '_gene.gff'

        s0_subfunctions.subfunction_make_acr_close_to_gene_cate(input_rice_gff_fl, input_rice_acr_fl,
                                                                s3_1_store_acr_cate_nearbyGene_dir, 'rice')
        s0_subfunctions.subfunction_make_acr_close_to_gene_cate(input_maize_gff_fl, input_maize_acr_fl,
                                                                s3_1_store_acr_cate_nearbyGene_dir, 'maize')
        s0_subfunctions.subfunction_make_acr_close_to_gene_cate(input_sorghum_gff_fl, input_sorghum_acr_fl,
                                                                s3_1_store_acr_cate_nearbyGene_dir, 'sorghum')

    if s3_s2_open_summarize_overview_cross_species == 'yes':

        s3_s2_open_summarize_overview_cross_species_dir = step03_summarize_overview_H3K27me3_dir + '/s3_s2_open_summarize_overview_cross_species_dir'
        if not os.path.exists(s3_s2_open_summarize_overview_cross_species_dir):
            os.makedirs(s3_s2_open_summarize_overview_cross_species_dir)

        s3_s2_target_spe2_list = s3_s2_target_spe2_str.split(',')

        for eachspe2 in s3_s2_target_spe2_list:

            ###############
            ##for the maize
            store_rice_to_maize_dir = s3_s2_open_summarize_overview_cross_species_dir + '/store_' + s3_s2_target_spe1 + '_to_' + eachspe2 + '_dir'
            if not os.path.exists(store_rice_to_maize_dir):
                os.makedirs(store_rice_to_maize_dir)

            spe1_prefix = s3_s2_target_spe1
            spe2_prefix = eachspe2

            ipt_spe1_acr_cate_fl = step03_summarize_overview_H3K27me3_dir + '/s3_1_store_acr_cate_nearbyGene_dir/opt_' + spe1_prefix + '_acr_toGeneCate.txt'
            ipt_spe2_acr_cate_fl = step03_summarize_overview_H3K27me3_dir + '/s3_1_store_acr_cate_nearbyGene_dir/opt_' + spe2_prefix + '_acr_toGeneCate.txt'

            ipt_final_summary_blast_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix \
                                         + '_dir/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'




            ipt_target_organ_H3K27_peak_rice_fl = input_target_organ_H3K27_peak_all_dir + '/' + spe1_prefix + '_H3K27me3_peak.txt'
            ipt_target_organ_H3K27_peak_maize_fl = input_target_organ_H3K27_peak_all_dir + '/' + spe2_prefix + '_H3K27me3_peak.txt'

            s3_subfunctions.subfunction_summarize_overview_for_H3K27me3(ipt_final_summary_blast_fl, ipt_target_organ_H3K27_peak_rice_fl,ipt_target_organ_H3K27_peak_maize_fl,
                                                                        ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,
                                                                        spe1_prefix, spe2_prefix, store_rice_to_maize_dir,
                                                                        s3_s2_extend_H3K27me3_flank_bp)

            ##updating 010324
            ipt_final_summary_blast_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix \
                                         + '_dir/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype_addSpeSpec.txt'

            store_rice_to_maize_dir = s3_s2_open_summarize_overview_cross_species_dir + '/store_' + s3_s2_target_spe1 + '_to_' + eachspe2 + '_include_species_specific_dir'
            if not os.path.exists(store_rice_to_maize_dir):
                os.makedirs(store_rice_to_maize_dir)

            ipt_target_organ_H3K27_peak_rice_fl = input_target_organ_H3K27_peak_all_dir + '/' + spe1_prefix + '_H3K27me3_peak.txt'
            ipt_target_organ_H3K27_peak_maize_fl = input_target_organ_H3K27_peak_all_dir + '/' + spe2_prefix + '_H3K27me3_peak.txt'

            s3_subfunctions.subfunction_summarize_overview_for_H3K27me3(ipt_final_summary_blast_fl,
                                                                        ipt_target_organ_H3K27_peak_rice_fl,
                                                                        ipt_target_organ_H3K27_peak_maize_fl,
                                                                        ipt_spe1_acr_cate_fl, ipt_spe2_acr_cate_fl,
                                                                        spe1_prefix, spe2_prefix,
                                                                        store_rice_to_maize_dir,
                                                                        s3_s2_extend_H3K27me3_flank_bp)


    if s3_s3_open_check_conserved_H3K27me3_broad_ACR_in_other_tissue == 'yes':

        s3_s3_target_organ_list = s3_s3_target_organ_str.split(',')

        s3_s3_open_check_conserved_H3K27me3_broad_ACR_in_other_tissue_dir = step03_summarize_overview_H3K27me3_dir + '/s3_s3_open_check_conserved_H3K27me3_broad_ACR_in_other_tissue_dir'
        if not os.path.exists(s3_s3_open_check_conserved_H3K27me3_broad_ACR_in_other_tissue_dir):
            os.makedirs(s3_s3_open_check_conserved_H3K27me3_broad_ACR_in_other_tissue_dir)

        s3_s2_target_spe2_list = s3_s2_target_spe2_str.split(',')

        for eachspe2 in s3_s2_target_spe2_list:

            store_rice_to_maize_dir = s3_s3_open_check_conserved_H3K27me3_broad_ACR_in_other_tissue_dir + '/store_' + s3_s2_target_spe1 + '_to_' + eachspe2 + '_dir'
            if not os.path.exists(store_rice_to_maize_dir):
                os.makedirs(store_rice_to_maize_dir)

            for eachorgan in s3_s3_target_organ_list:

                spe1_prefix = s3_s2_target_spe1
                spe2_prefix = eachspe2

                ipt_other_organ_fl = input_rice_all_organ_H3K27_fl_dir + '/opt_rice_' + eachorgan + '_H3K27me3.txt'

                ipt_final_summary_blast_fl = step03_summarize_overview_H3K27me3_dir + '/s3_s2_open_summarize_overview_cross_species_dir/' + 'store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/' + 'opt_final_summary_file_add_' + spe1_prefix + '_' + spe2_prefix + '.txt'

                s3_subfunctions.subfunction_check_conserved_H3K27me3_in_other_organs (ipt_final_summary_blast_fl,ipt_other_organ_fl,
                                                                                      store_rice_to_maize_dir,spe1_prefix,spe2_prefix,eachorgan)


    if s3_s4_open_check_motif_freq_in_diff_cate_acr == 'yes':

        s3_s4_open_check_motif_freq_in_diff_cate_acr_dir = step03_summarize_overview_H3K27me3_dir + '/s3_s4_open_check_motif_freq_in_diff_cate_acr_dir'
        if not os.path.exists(s3_s4_open_check_motif_freq_in_diff_cate_acr_dir):
            os.makedirs(s3_s4_open_check_motif_freq_in_diff_cate_acr_dir)

        s3_s2_target_spe2_list = s3_s2_target_spe2_str.split(',')

        for eachspe2 in s3_s2_target_spe2_list:

            store_rice_to_maize_dir = s3_s4_open_check_motif_freq_in_diff_cate_acr_dir + '/store_' + s3_s2_target_spe1 + '_to_' + eachspe2 + '_dir'
            if not os.path.exists(store_rice_to_maize_dir):
                os.makedirs(store_rice_to_maize_dir)

            spe1_prefix = s3_s2_target_spe1
            spe2_prefix = eachspe2

            ipt_final_summary_fl = step03_summarize_overview_H3K27me3_dir + '/s3_s2_open_summarize_overview_cross_species_dir/' + 'store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/' + 'opt_final_summary_file_add_' + spe1_prefix + '_' + spe2_prefix + '.txt'

            ipt_motif_fimo_fl = input_all_spe_PRE_motif_fimo_dir + '/opt_' + s3_s2_target_spe1 + '_PRE_motif_fimo.txt'
            opt_dir = store_rice_to_maize_dir
            s3_subfunctions.subfunction_check_motif_freq_in_diff_cate(ipt_final_summary_fl, ipt_motif_fimo_fl, opt_dir, spe1_prefix,
                                                      spe2_prefix)


    if s3_s5_open_check_cns_overlap_ACR == 'yes':

        s3_s5_open_check_cns_overlap_ACR_dir = step03_summarize_overview_H3K27me3_dir + '/s3_s5_open_check_cns_overlap_ACR_dir'
        if not os.path.exists(s3_s5_open_check_cns_overlap_ACR_dir):
            os.makedirs(s3_s5_open_check_cns_overlap_ACR_dir)

        s3_s2_target_spe2_list = s3_s2_target_spe2_str.split(',')

        for eachspe2 in s3_s2_target_spe2_list:

            store_rice_to_maize_dir = s3_s5_open_check_cns_overlap_ACR_dir + '/store_' + s3_s2_target_spe1 + '_to_' + eachspe2 + '_dir'
            if not os.path.exists(store_rice_to_maize_dir):
                os.makedirs(store_rice_to_maize_dir)

            spe1_prefix = s3_s2_target_spe1
            spe2_prefix = eachspe2

            #ipt_final_summary_fl = step03_summarize_overview_H3K27me3_dir + '/s3_s2_open_summarize_overview_cross_species_dir/' + 'store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/' + 'opt_final_summary_file_add_' + spe1_prefix + '_' + spe2_prefix + '.txt'

            #ipt_spe1_cns_gff_fl = input_all_spe_cns_gff_fl_dir + '/ipt_' + spe1_prefix + '_cns_gff.txt'
            ##opt_dir = store_rice_to_maize_dir
            #s3_subfunctions.subfunction_overlap_with_CNS (ipt_final_summary_fl,ipt_spe1_cns_gff_fl,spe1_prefix,spe2_prefix,opt_dir)


            ##updating 010824
            ##we will include final spe specific to have a check
            ipt_final_summary_fl = step03_summarize_overview_H3K27me3_dir + '/s3_s2_open_summarize_overview_cross_species_dir' + '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_include_species_specific_dir/' + \
                                   'opt_final_summary_file_add_' + spe1_prefix + '_' + spe2_prefix + '.txt'
            ipt_spe1_cns_gff_fl = input_all_spe_cns_gff_fl_dir + '/ipt_' + spe1_prefix + '_cns_gff.txt'
            opt_dir = store_rice_to_maize_dir
            s3_subfunctions.subfunction_overlap_with_CNS(ipt_final_summary_fl, ipt_spe1_cns_gff_fl, spe1_prefix,
                                                         spe2_prefix, opt_dir)


    if s3_s6_open_check_overlap_atlas_acr == 'yes':

        s3_s6_open_check_overlap_atlas_acr = step03_summarize_overview_H3K27me3_dir + '/s3_s6_open_check_overlap_atlas_acr'
        if not os.path.exists(s3_s6_open_check_overlap_atlas_acr):
            os.makedirs(s3_s6_open_check_overlap_atlas_acr)

        s3_s2_target_spe2_list = s3_s2_target_spe2_str.split(',')

        for eachspe2 in s3_s2_target_spe2_list:

            store_rice_to_maize_dir = s3_s6_open_check_overlap_atlas_acr + '/store_' + s3_s2_target_spe1 + '_to_' + eachspe2 + '_dir'
            if not os.path.exists(store_rice_to_maize_dir):
                os.makedirs(store_rice_to_maize_dir)

            spe1_prefix = s3_s2_target_spe1
            spe2_prefix = eachspe2

            ipt_rice_atlas_fl = input_total_rice_atlas_acr_fl
            ipt_final_summary_fl = step03_summarize_overview_H3K27me3_dir + '/s3_s2_open_summarize_overview_cross_species_dir/' + 'store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/' + 'opt_final_summary_file_add_' + spe1_prefix + '_' + spe2_prefix + '.txt'

            opt_dir = store_rice_to_maize_dir

            s3_subfunctions.subfunction_overlap_with_atlas_acr (ipt_rice_atlas_fl,input_rice_atlas_acr_coverage_fl,ipt_final_summary_fl,
                                                                input_rice_all_organ_H3K27_fl_dir,
                                                                s2_s2_organct_coverage_cutoff,opt_dir,spe1_prefix,spe2_prefix)


    if s3_s7_open_check_target_motif_in_ACR == 'yes':

        s3_s7_open_check_target_motif_in_ACR_dir = step03_summarize_overview_H3K27me3_dir + '/s3_s7_open_check_target_motif_in_ACR_dir'
        if not os.path.exists(s3_s7_open_check_target_motif_in_ACR_dir):
            os.makedirs(s3_s7_open_check_target_motif_in_ACR_dir)

        s3_s2_target_spe2_list = s3_s2_target_spe2_str.split(',')

        for eachspe2 in s3_s2_target_spe2_list:

            store_rice_to_maize_dir = s3_s7_open_check_target_motif_in_ACR_dir + '/store_' + s3_s2_target_spe1 + '_to_' + eachspe2 + '_dir'
            if not os.path.exists(store_rice_to_maize_dir):
                os.makedirs(store_rice_to_maize_dir)

            spe1_prefix = s3_s2_target_spe1
            spe2_prefix = eachspe2

            ipt_final_summary_fl = step03_summarize_overview_H3K27me3_dir + '/s3_s2_open_summarize_overview_cross_species_dir' + \
                                   '/store_' + spe1_prefix + '_to_' + spe2_prefix + '_include_species_specific_dir/' + \
                                   'opt_final_summary_file_add_' + spe1_prefix + '_' + spe2_prefix + '.txt'

            ipt_spe1_motif_fl = input_all_motif_fimo_dir + '/' + spe1_prefix + '_motif_fimo.txt'
            ipt_spe2_motif_fl = input_all_motif_fimo_dir + '/' + spe2_prefix + '_motif_fimo.txt'

            s3_subfunctions.subfunction_add_target_motif_in_final_summary(ipt_final_summary_fl, s3_s7_target_motifID,
                                                                          ipt_spe1_motif_fl,ipt_spe2_motif_fl,store_rice_to_maize_dir,spe1_prefix, spe2_prefix)






########
##step04
def step04_intersect_three_species (input_output_dir,s0_s1_target_spe1,
                                    s4_s1_open_find_candidate_plot_acr_H3K27me3,
                                    s4_s2_open_find_candidate_plot_five_spe_acr,
                                    s4_s3_open_find_candidate_plot_target_celltype,s4_s3_target_spe2_str,s4_s3_target_celltype_str):

    step04_intersect_three_species_dir = input_output_dir + '/step04_intersect_three_species_dir'
    if not os.path.exists(step04_intersect_three_species_dir):
        os.makedirs(step04_intersect_three_species_dir)

    if s4_s1_open_find_candidate_plot_acr_H3K27me3 == 'yes':

        ##subfunction_iden_shared_acr_under_H3K27me3 (ipt_spe1_to_spe2_final_fl,ipt_spe1_to_spe3_final_fl,ipt_spe1_acr_cate_fl,ipt_spe2_acr_cate_fl,ipt_spe3_acr_cate_fl,opt_dir)
        s4_s1_open_find_candidate_plot_acr_H3K27me3_dir = step04_intersect_three_species_dir + '/s4_s1_open_find_candidate_plot_acr_H3K27me3_dir'
        if not os.path.exists(s4_s1_open_find_candidate_plot_acr_H3K27me3):
            os.makedirs(s4_s1_open_find_candidate_plot_acr_H3K27me3)

        spe1_prefix = 'rice'
        spe2_prefix = 'maize'
        ipt_spe1_to_spe2_final_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

        spe1_prefix = 'rice'
        spe3_prefix = 'sorghum'
        ipt_spe1_to_spe3_final_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe3_prefix + '_dir/opt_' + spe1_prefix + '_' + spe3_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

        #'/opt_' + ipt_prefix + '_acr_toGeneCate.txt'
        ipt_spe1_acr_cate_fl = input_output_dir + '/step02_summarize_overview_dir/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/opt_' + spe1_prefix + '_acr_toGeneCate.txt'
        ipt_spe2_acr_cate_fl = input_output_dir + '/step02_summarize_overview_dir/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/opt_' + spe2_prefix + '_acr_toGeneCate.txt'
        ipt_spe3_acr_cate_fl = input_output_dir + '/step02_summarize_overview_dir/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/opt_' + spe3_prefix + '_acr_toGeneCate.txt'

        store_final_overlapped_acr_line_BroadToBroad_H3K27me3_list, store_final_overlapped_acr_line_BroadToRestricted_H3K27me3_list = \
            s4_subfunctions.subfunction_iden_shared_acr_under_H3K27me3(ipt_spe1_to_spe2_final_fl,ipt_spe1_to_spe3_final_fl,ipt_spe1_acr_cate_fl,
                                                                   ipt_spe2_acr_cate_fl,ipt_spe3_acr_cate_fl,s4_s1_open_find_candidate_plot_acr_H3K27me3_dir,spe1_prefix,spe2_prefix,spe3_prefix)

        with open(s4_s1_open_find_candidate_plot_acr_H3K27me3 + '/opt_final_intersect_candidate_acr_BroadToBroad_H3K27me3_plot.txt', 'w+') as opt:
            for eachline in store_final_overlapped_acr_line_BroadToBroad_H3K27me3_list:
                opt.write(eachline + '\n')

        with open(s4_s1_open_find_candidate_plot_acr_H3K27me3 + '/opt_final_intersect_candidate_acr_BroadToRestricted_H3K27me3_plot.txt', 'w+') as opt:
            for eachline in store_final_overlapped_acr_line_BroadToRestricted_H3K27me3_list:
                opt.write(eachline + '\n')


    if s4_s2_open_find_candidate_plot_five_spe_acr == 'yes':

        s4_s2_open_find_candidate_plot_five_spe_acr_dir = step04_intersect_three_species_dir + '/s4_s2_open_find_candidate_plot_five_spe_acr_dir'
        if not os.path.exists(s4_s2_open_find_candidate_plot_five_spe_acr_dir):
            os.makedirs(s4_s2_open_find_candidate_plot_five_spe_acr_dir)

        spe1_prefix = 'rice'
        spe2_prefix = 'maize'
        ipt_spe1_to_spe2_final_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

        spe1_prefix = 'rice'
        spe3_prefix = 'sorghum'
        ipt_spe1_to_spe3_final_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe3_prefix + '_dir/opt_' + spe1_prefix + '_' + spe3_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

        spe1_prefix = 'rice'
        spe4_prefix = 'Pm'
        ipt_spe1_to_spe4_final_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe4_prefix + '_dir/opt_' + spe1_prefix + '_' + spe4_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

        spe1_prefix = 'rice'
        spe5_prefix = 'Uf'
        ipt_spe1_to_spe5_final_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + spe1_prefix + '_to_' + spe5_prefix + '_dir/opt_' + spe1_prefix + '_' + spe5_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

        ipt_spe1_acr_cate_fl = input_output_dir + '/step02_summarize_overview_dir/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/opt_' + spe1_prefix + '_acr_toGeneCate.txt'
        ipt_spe2_acr_cate_fl = input_output_dir + '/step02_summarize_overview_dir/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/opt_' + spe2_prefix + '_acr_toGeneCate.txt'
        ipt_spe3_acr_cate_fl = input_output_dir + '/step02_summarize_overview_dir/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/opt_' + spe3_prefix + '_acr_toGeneCate.txt'
        ipt_spe4_acr_cate_fl = input_output_dir + '/step02_summarize_overview_dir/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/opt_' + spe4_prefix + '_acr_toGeneCate.txt'
        ipt_spe5_acr_cate_fl = input_output_dir + '/step02_summarize_overview_dir/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/opt_' + spe5_prefix + '_acr_toGeneCate.txt'

        store_final_overlapped_acr_line_five_species_list = s4_subfunctions.subfunction_iden_shared_acr_five_species(ipt_spe1_to_spe2_final_fl, ipt_spe1_to_spe3_final_fl,
                                                 ipt_spe1_to_spe4_final_fl, ipt_spe1_to_spe5_final_fl,
                                                 ipt_spe1_acr_cate_fl, ipt_spe2_acr_cate_fl, ipt_spe3_acr_cate_fl,
                                                 ipt_spe4_acr_cate_fl, ipt_spe5_acr_cate_fl,
                                                 s4_s2_open_find_candidate_plot_five_spe_acr_dir,
                                                 spe1_prefix, spe2_prefix, spe3_prefix, spe4_prefix,
                                                 spe5_prefix)



        with open(s4_s2_open_find_candidate_plot_five_spe_acr_dir + '/opt_final_intersect_candidate_acr_H3K27me3_plot.txt', 'w+') as opt:
            for eachline in store_final_overlapped_acr_line_five_species_list:
                opt.write(eachline + '\n')

    if s4_s3_open_find_candidate_plot_target_celltype == 'yes':

        s4_s3_open_find_candidate_plot_target_celltype_dir = step04_intersect_three_species_dir + '/s4_s3_open_find_candidate_plot_target_celltype_dir'
        if not os.path.exists(s4_s3_open_find_candidate_plot_target_celltype_dir):
            os.makedirs(s4_s3_open_find_candidate_plot_target_celltype_dir)

        s4_s3_target_spe2_list = s4_s3_target_spe2_str.split(',')
        s4_s3_target_celltype_list = s4_s3_target_celltype_str.split(',')

        for eachspe2 in s4_s3_target_spe2_list:

            opt_target_spe_dir = s4_s3_open_find_candidate_plot_target_celltype_dir + '/' + eachspe2
            if not os.path.exists(opt_target_spe_dir):
                os.makedirs(opt_target_spe_dir)

            for eachcelltype in s4_s3_target_celltype_list:

                spe1_prefix = s0_s1_target_spe1
                spe2_prefix = eachspe2

                ipt_spe1_to_spe2_final_fl = input_output_dir + '/step01_species_compare_add_cate_dir/store_' + \
                                            spe1_prefix + '_to_' + spe2_prefix + '_dir/opt_' + spe1_prefix + '_' + spe2_prefix + '_H3K27me3_addCelltype_SynRegion_addDirFlt_addRegFlt_addspe2DirFlt_addspe2CT_addCoverCelltype.txt'

                ipt_spe1_acr_cate_fl = input_output_dir + '/step02_summarize_overview_dir/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/opt_' + spe1_prefix + '_acr_toGeneCate.txt'
                ipt_spe2_acr_cate_fl = input_output_dir + '/step02_summarize_overview_dir/s2_s1_open_basic_summary_overview_dir/store_acr_cate_nearbyGene_dir/opt_' + spe2_prefix + '_acr_toGeneCate.txt'

                store_final_overlapped_acr_line_celltype_list = s4_subfunctions.subfunction_iden_celltype_specific_ACR_correspond(ipt_spe1_to_spe2_final_fl,
                                                                  ipt_spe1_acr_cate_fl, ipt_spe2_acr_cate_fl,
                                                                  spe1_prefix, spe2_prefix,
                                                                  eachcelltype, opt_target_spe_dir)

                with open(opt_target_spe_dir + '/opt_final_intersect_candidate_acr_' + eachcelltype + '_plot.txt',
                        'w+') as opt:
                    for eachline in store_final_overlapped_acr_line_celltype_list:
                        opt.write(eachline + '\n')








store_target_parameter_dic = {}
with open (input_configure_fl,'r') as ipt:
    for eachline in ipt:
        eachline = eachline.strip('\n')
        if not eachline.startswith('#'):
            col = eachline.strip().split('=')
            store_target_parameter_dic[col[0]] = col[1]

step00 = store_target_parameter_dic['step00'] ##yes or no
step01 = store_target_parameter_dic['step01'] ##yes or no
step02 = store_target_parameter_dic['step02'] ##yes or no
step03 = store_target_parameter_dic['step03'] ##yes or no
step04 = store_target_parameter_dic['step04']

input_core_num = store_target_parameter_dic['input_core_num']

s0_s1_open_check_syntenic_block = store_target_parameter_dic['s0_s1_open_check_syntenic_block']
s0_s1_target_spe1 = store_target_parameter_dic['s0_s1_target_spe1']
s0_s1_target_spe2_str = store_target_parameter_dic['s0_s1_target_spe2_str']
s0_s1_width_bin = store_target_parameter_dic['s0_s1_width_bin']

s0_s2_open_check_te_composition = store_target_parameter_dic['s0_s2_open_check_te_composition']
s0_s3_open_check_gene_cate = store_target_parameter_dic['s0_s3_open_check_gene_cate']

s0_s4_open_motif_coverage = store_target_parameter_dic['s0_s4_open_motif_coverage']
s0_s4_target_species_str = store_target_parameter_dic['s0_s4_target_species_str']
s0_s4_sub_open_build_control = store_target_parameter_dic['s0_s4_sub_open_build_control']
s0_s4_sub_open_check_avg_for_acr = store_target_parameter_dic['s0_s4_sub_open_check_avg_for_acr']
s0_s4_sub_extend_acr_range_bp = store_target_parameter_dic['s0_s4_sub_extend_acr_range_bp']
s0_s4_sub_open_build_posi_coverage = store_target_parameter_dic['s0_s4_sub_open_build_posi_coverage']

s0_s5_open_motifs_TFs = store_target_parameter_dic['s0_s5_open_motifs_TFs']
s0_s5_sub_open_prepare_gene_accessibility = store_target_parameter_dic['s0_s5_sub_open_prepare_gene_accessibility']
s0_s5_sub_target_spe_str_for_geneacc_smooth = store_target_parameter_dic['s0_s5_sub_target_spe_str_for_geneacc_smooth']
s0_s5_sub_open_smooth_gene = store_target_parameter_dic['s0_s5_sub_open_smooth_gene']
s0_s5_sub_target_cluster = store_target_parameter_dic['s0_s5_sub_target_cluster']

s0_s5_sub_open_smooth_motif = store_target_parameter_dic['s0_s5_sub_open_smooth_motif']
s0_s5_sub_open_corr_motif_TF = store_target_parameter_dic['s0_s5_sub_open_corr_motif_TF']

s0_s5_sub_open_plot_TF_motif = store_target_parameter_dic['s0_s5_sub_open_plot_TF_motif']
s0_s5_sub_lim = store_target_parameter_dic['s0_s5_sub_lim']
s0_s5_sub_plot_topnum_pairs = store_target_parameter_dic['s0_s5_sub_plot_topnum_pairs']

s0_s5_sub_open_plot_TFfam_TF_motif = store_target_parameter_dic['s0_s5_sub_open_plot_TFfam_TF_motif']
s0_s5_sub_target_fam_str = store_target_parameter_dic['s0_s5_sub_target_fam_str']

s0_s5_sub_open_plot_high_corr_all_spe = store_target_parameter_dic['s0_s5_sub_open_plot_high_corr_all_spe']
s0_s5_sub_target_plot_spe_str = store_target_parameter_dic['s0_s5_sub_target_plot_spe_str']
s0_s5_sub_target_plot_corr_cutoff = store_target_parameter_dic['s0_s5_sub_target_plot_corr_cutoff']

##updating 123023
s0_s5_sub_open_plot_Uf_TF_motif = store_target_parameter_dic['s0_s5_sub_open_plot_Uf_TF_motif']
s0_s5_sub_motif_TF_str_for_Uf = store_target_parameter_dic['s0_s5_sub_motif_TF_str_for_Uf']

s0_s6_model_based_motif_enrich = store_target_parameter_dic['s0_s6_model_based_motif_enrich']
s0_s6_sub_open_flt_meta_cell = store_target_parameter_dic['s0_s6_sub_open_flt_meta_cell']
s0_s6_sub_cluster_colnum_str = store_target_parameter_dic['s0_s6_sub_cluster_colnum_str']

s0_s6_sub_target_spe_str_for_add_ACRnum = store_target_parameter_dic['s0_s6_sub_target_spe_str_for_add_ACRnum']
s0_s6_sub_open_add_ACRnum_to_meta = store_target_parameter_dic['s0_s6_sub_open_add_ACRnum_to_meta']

s0_s6_sub_open_add_log10tn5_lib_to_meta = store_target_parameter_dic['s0_s6_sub_open_add_log10tn5_lib_to_meta']
s0_s6_sub_target_library_colnum = store_target_parameter_dic['s0_s6_sub_target_library_colnum'] ##18
s0_s6_sub_target_tn5_colnum = store_target_parameter_dic['s0_s6_sub_target_tn5_colnum'] ##2

s0_s6_sub_open_prepare_acrmotif = store_target_parameter_dic['s0_s6_sub_open_prepare_acrmotif']
s0_s6_sub_target_spe_str_for_conduct_enrich_str = store_target_parameter_dic['s0_s6_sub_target_spe_str_for_conduct_enrich_str']
s0_s6_sub_open_prepare_enrich = store_target_parameter_dic['s0_s6_sub_open_prepare_enrich']
s0_s6_sub_pvalcutoff = store_target_parameter_dic['s0_s6_sub_pvalcutoff']
s0_s6_sub_open_conduct_enrich = store_target_parameter_dic['s0_s6_sub_open_conduct_enrich']
s0_s6_sub_target_clusternm = store_target_parameter_dic['s0_s6_sub_target_clusternm']

##updating 121823
s0_s7_open_check_cns = store_target_parameter_dic['s0_s7_open_check_cns']
s0_s7_sub_open_eachspe_overlapCNS = store_target_parameter_dic['s0_s7_sub_open_eachspe_overlapCNS']
s0_s7_sub_target_all_spe_str = store_target_parameter_dic['s0_s7_sub_target_all_spe_str']
s0_s7_sub_open_rice_atlas_CNS = store_target_parameter_dic['s0_s7_sub_open_rice_atlas_CNS']


s2_s1_open_basic_summary_overview = store_target_parameter_dic['s2_s1_open_basic_summary_overview']
s2_s1_open_iden_acr_to_gene_cate = store_target_parameter_dic['s2_s1_open_iden_acr_to_gene_cate']
s2_s1_open_check_phy_score = store_target_parameter_dic['s2_s1_open_check_phy_score']

s2_s2_oepn_check_atlasACR_celltypenum = store_target_parameter_dic['s2_s2_oepn_check_atlasACR_celltypenum']
s2_s2_organct_coverage_cutoff = store_target_parameter_dic['s2_s2_organct_coverage_cutoff']

s2_s3_sub_open_only_add_spe1_cpm = store_target_parameter_dic['s2_s3_sub_open_only_add_spe1_cpm']
s2_s3_open_check_nearby_gene = store_target_parameter_dic['s2_s3_open_check_nearby_gene']
s2_s3_target_celltype_addCpm_str = store_target_parameter_dic['s2_s3_target_celltype_addCpm_str']

s2_s3_sub_open_add_spe1_spe2_cpm = store_target_parameter_dic['s2_s3_sub_open_add_spe1_spe2_cpm']
s2_s3_target_spe1_for_twospe_cpm = store_target_parameter_dic['s2_s3_target_spe1_for_twospe_cpm']
s2_s3_target_spe2_str_for_twospe_cpm = store_target_parameter_dic['s2_s3_target_spe2_str_for_twospe_cpm']
s2_s3_target_celltype_addCpm_spe2_str = store_target_parameter_dic['s2_s3_target_celltype_addCpm_spe2_str']

s2_s3_sub_open_add_celltype_specific_str = store_target_parameter_dic['s2_s3_sub_open_add_celltype_specific_str']
s2_s3_target_spe1_for_add_ct_str = store_target_parameter_dic['s2_s3_target_spe1_for_add_ct_str']
s2_s3_target_spe2_str_for_add_ct_str = store_target_parameter_dic['s2_s3_target_spe2_str_for_add_ct_str']
s2_s3_sub_cutoff_avg_fold = store_target_parameter_dic['s2_s3_sub_cutoff_avg_fold']


s2_s4_open_conduct_enrich_cate_celltype = store_target_parameter_dic['s2_s4_open_conduct_enrich_cate_celltype']
s2_s4_open_enrich = store_target_parameter_dic['s2_s4_open_enrich']
s2_s4_repeat_times = store_target_parameter_dic['s2_s4_repeat_times']

s2_s5_open_identify_target_motif_binding_plot = store_target_parameter_dic['s2_s5_open_identify_target_motif_binding_plot']
s2_s5_target_motif_id_str = store_target_parameter_dic['s2_s5_target_motif_id_str']


s2_s6_open_check_cns_num_per_acr = store_target_parameter_dic['s2_s6_open_check_cns_num_per_acr']
s2_s6_target_spe1 = store_target_parameter_dic['s2_s6_target_spe1']
s2_s6_target_spe2_str = store_target_parameter_dic['s2_s6_target_spe2_str']

s2_s6_sub_open_check_CNS_overlap_ACR = store_target_parameter_dic['s2_s6_sub_open_check_CNS_overlap_ACR']
s2_s6_sub_open_varACRs_in_otherSpes_blastToregionInrice_shown_in_atlas = store_target_parameter_dic['s2_s6_sub_open_varACRs_in_otherSpes_blastToregionInrice_shown_in_atlas']

s2_s7_open_check_snp_num_per_acr = store_target_parameter_dic['s2_s7_open_check_snp_num_per_acr']
s2_s7_target_spe1 = store_target_parameter_dic['s2_s7_target_spe1']
s2_s7_target_spe2_str = store_target_parameter_dic['s2_s7_target_spe2_str']

s2_s8_open_motif_enrich_nonsyn = store_target_parameter_dic['s2_s8_open_motif_enrich_nonsyn']
s2_s8_target_spe1 = store_target_parameter_dic['s2_s8_target_spe1']
s2_s8_target_spe2_str = store_target_parameter_dic['s2_s8_target_spe2_str']

s2_s8_sub_method = store_target_parameter_dic['s2_s8_sub_method']

s2_s9_open_check_H3K27m3_under_species_specific_ACR = store_target_parameter_dic['s2_s9_open_check_H3K27m3_under_species_specific_ACR']
s2_s9_target_spe1 = store_target_parameter_dic['s2_s9_target_spe1']
s2_s9_target_spe2_str = store_target_parameter_dic['s2_s9_target_spe2_str']

s2_s10_open_check_motif_enrich_on_exp_cor_ACR = store_target_parameter_dic['s2_s10_open_check_motif_enrich_on_exp_cor_ACR']
s2_s10_target_spe1 = store_target_parameter_dic['s2_s10_target_spe1']
s2_s10_target_ct_str = store_target_parameter_dic['s2_s10_target_ct_str']
s2_s10_repeat_times = store_target_parameter_dic['s2_s10_repeat_times']

s2_s11_open_check_TE_on_ACR = store_target_parameter_dic['s2_s11_open_check_TE_on_ACR']

s2_s11_sub_open_overlap_enrich_TE = store_target_parameter_dic['s2_s11_sub_open_overlap_enrich_TE']
s2_s11_target_spe1 = store_target_parameter_dic['s2_s11_target_spe1']
s2_s11_target_spe2_str = store_target_parameter_dic['s2_s11_target_spe2_str']
s2_s11_sub_open_motif_enrich_target_TEfam = store_target_parameter_dic['s2_s11_sub_open_motif_enrich_target_TEfam']
s2_s11_sub_target_CT_TEfam_str = store_target_parameter_dic['s2_s11_sub_target_CT_TEfam_str']
s2_s11_sub_repeat_times = store_target_parameter_dic['s2_s11_sub_repeat_times']

s2_s11_sub_open_DNAmethylation_target_TEfam = store_target_parameter_dic['s2_s11_sub_open_DNAmethylation_target_TEfam']


s3_s1_open_iden_acr_to_gene_cate = store_target_parameter_dic['s3_s1_open_iden_acr_to_gene_cate']
s3_s1_target_spe1 = store_target_parameter_dic['s3_s1_target_spe1']
s3_s1_target_spe2_str = store_target_parameter_dic['s3_s1_target_spe2_str']
s3_s2_open_summarize_overview_cross_species = store_target_parameter_dic['s3_s2_open_summarize_overview_cross_species']
s3_s2_extend_H3K27me3_flank_bp = store_target_parameter_dic['s3_s2_extend_H3K27me3_flank_bp']
s3_s2_target_spe1 = store_target_parameter_dic['s3_s2_target_spe1']
s3_s2_target_spe2_str = store_target_parameter_dic['s3_s2_target_spe2_str']
s3_s3_open_check_conserved_H3K27me3_broad_ACR_in_other_tissue = store_target_parameter_dic['s3_s3_open_check_conserved_H3K27me3_broad_ACR_in_other_tissue']
s3_s3_target_organ_str = store_target_parameter_dic['s3_s3_target_organ_str']

#s3_open_iden_acr_to_gene_cate = store_target_parameter_dic['s3_open_iden_acr_to_gene_cate']
#s3_extend_H3K27me3_flank_bp = store_target_parameter_dic['s3_extend_H3K27me3_flank_bp']


s3_s4_open_check_motif_freq_in_diff_cate_acr = store_target_parameter_dic['s3_s4_open_check_motif_freq_in_diff_cate_acr']

s3_s5_open_check_cns_overlap_ACR = store_target_parameter_dic['s3_s5_open_check_cns_overlap_ACR']

s3_s6_open_check_overlap_atlas_acr = store_target_parameter_dic['s3_s6_open_check_overlap_atlas_acr']


s3_s7_open_check_target_motif_in_ACR = store_target_parameter_dic['s3_s7_open_check_target_motif_in_ACR']
s3_s7_target_motifID = store_target_parameter_dic['s3_s7_target_motifID']


s4_s1_open_find_candidate_plot_acr_H3K27me3 = store_target_parameter_dic['s4_s1_open_find_candidate_plot_acr_H3K27me3']
s4_s2_open_find_candidate_plot_five_spe_acr = store_target_parameter_dic['s4_s2_open_find_candidate_plot_five_spe_acr']

#if step001 == 'yes':
s4_s3_open_find_candidate_plot_target_celltype = store_target_parameter_dic['s4_s3_open_find_candidate_plot_target_celltype']
s4_s3_target_spe2_str = store_target_parameter_dic['s4_s3_target_spe2_str']
s4_s3_target_celltype_str = store_target_parameter_dic['s4_s3_target_celltype_str']



if step00 == 'yes':
    step00_basic_characters(input_required_scripts_dir, input_output_dir, input_rice_syntenic_genes_all_os_ACRs_fl,
                            input_all_syntenic_regions_os_acrs_dir,
                            input_all_H3K27me3_dir, input_rice_acr_add_celltype_fl, input_total_rice_atlas_acr_fl,
                            input_all_spe_TE_dir, input_all_celltype_acr_dir, input_all_spe_gene_gff_dir,
                            input_all_genome_size_dir,
                            input_all_motif_fimo_dir,
                            input_gene_raw_sparse_fl_dir, input_all_spe_meta_fl_dir, input_all_spe_SVD_fl_dir,
                            input_all_spe_motif_dev_score_dir,
                            input_all_spe_ortho_TF_to_AtTF_dir, input_At_TF_to_motif_fl, input_Atgene_TFid_fl,
                            input_all_spe_peak_sparse_dir, ipt_motif_commonnm_fl,
                            input_all_spe_cns_gff_fl_dir,
                            input_core_num,
                            s0_s1_open_check_syntenic_block, s0_s1_target_spe1, s0_s1_target_spe2_str, s0_s1_width_bin,
                            s0_s2_open_check_te_composition,
                            s0_s3_open_check_gene_cate,
                            s0_s4_open_motif_coverage, s0_s4_target_species_str, s0_s4_sub_open_build_control,
                            s0_s4_sub_open_check_avg_for_acr, s0_s4_sub_extend_acr_range_bp,
                            s0_s4_sub_open_build_posi_coverage,
                            s0_s5_open_motifs_TFs,
                            s0_s5_sub_open_prepare_gene_accessibility, s0_s5_sub_target_spe_str_for_geneacc_smooth,
                            s0_s5_sub_open_smooth_gene, s0_s5_sub_target_cluster,
                            s0_s5_sub_open_smooth_motif,
                            s0_s5_sub_open_corr_motif_TF,
                            s0_s5_sub_open_plot_TF_motif, s0_s5_sub_lim, s0_s5_sub_plot_topnum_pairs,
                            s0_s5_sub_open_plot_TFfam_TF_motif, s0_s5_sub_target_fam_str,
                            s0_s5_sub_open_plot_high_corr_all_spe, s0_s5_sub_target_plot_spe_str,
                            s0_s5_sub_target_plot_corr_cutoff,
                            s0_s5_sub_open_plot_Uf_TF_motif, s0_s5_sub_motif_TF_str_for_Uf,
                            s0_s6_model_based_motif_enrich, s0_s6_sub_open_flt_meta_cell, s0_s6_sub_cluster_colnum_str,
                            s0_s6_sub_target_spe_str_for_add_ACRnum,
                            s0_s6_sub_open_add_ACRnum_to_meta,
                            s0_s6_sub_open_add_log10tn5_lib_to_meta, s0_s6_sub_target_library_colnum,
                            s0_s6_sub_target_tn5_colnum,
                            s0_s6_sub_open_prepare_acrmotif, s0_s6_sub_target_spe_str_for_conduct_enrich_str,
                            s0_s6_sub_open_prepare_enrich, s0_s6_sub_pvalcutoff,
                            s0_s6_sub_open_conduct_enrich, s0_s6_sub_target_clusternm,
                            s0_s7_open_check_cns, s0_s7_sub_open_eachspe_overlapCNS, s0_s7_sub_target_all_spe_str,
                            s0_s7_sub_open_rice_atlas_CNS)

if step01 == 'yes':
    step01_species_compare_add_cate(input_rice_to_allspe_blast_dir,
                                    input_all_celltype_acr_dir,
                                    input_all_spe_gene_gff_dir,
                                    input_rice_to_allspe_blast_record_riceID_dir,
                                    input_rice_syntenic_genes_all_os_ACRs_fl,
                                    input_all_syntenic_genes_all_os_ACRs_dir,
                                    input_all_H3K27me3_dir,
                                    input_output_dir, s0_s1_target_spe1, s0_s1_target_spe2_str)
#ipt_cpm_snRNAseq_fl ipt_maize_cpm_snRNAseq_fl
if step02 == 'yes':
    step02_summarize_overview(input_output_dir, input_all_syntenic_regions_os_acrs_dir, input_rice_phy_score_fl,
                              input_all_spe_gene_gff_dir, input_all_celltype_acr_dir, input_rice_atlas_acr_coverage_fl,
                              input_total_rice_atlas_acr_fl,
                              input_all_spe_gene_annot_fl_dir, input_all_motif_fimo_dir, input_all_spe_cns_gff_fl_dir,
                              input_marker_gene_list_fl_dir,
                              ipt_cpm_snRNAseq_fl, ipt_maize_cpm_snRNAseq_fl, ipt_cluster_annot_fl,
                              input_all_spe_geneID_symbol_dir, input_all_syntenic_genes_all_os_ACRs_dir,
                              input_all_spe_ortho_dir, input_H3K27me3_fl,
                              input_all_spe_TE_dir, input_DNAmethylation_bw_dir,
                              input_core_num,
                              s0_s1_target_spe1, s0_s1_target_spe2_str,
                              s2_s1_open_basic_summary_overview,
                              s2_s1_open_iden_acr_to_gene_cate,
                              s2_s1_open_check_phy_score,
                              s2_s2_oepn_check_atlasACR_celltypenum,
                              s2_s2_organct_coverage_cutoff,
                              s2_s3_sub_open_only_add_spe1_cpm, s2_s3_open_check_nearby_gene,
                              s2_s3_target_celltype_addCpm_str,
                              s2_s3_sub_open_add_spe1_spe2_cpm, s2_s3_target_spe1_for_twospe_cpm,
                              s2_s3_target_spe2_str_for_twospe_cpm,
                              s2_s3_target_celltype_addCpm_spe2_str,
                              s2_s3_sub_open_add_celltype_specific_str, s2_s3_target_spe1_for_add_ct_str,
                              s2_s3_target_spe2_str_for_add_ct_str, s2_s3_sub_cutoff_avg_fold,
                              s2_s4_open_conduct_enrich_cate_celltype, s2_s4_repeat_times, s2_s4_open_enrich,
                              s2_s5_open_identify_target_motif_binding_plot, s2_s5_target_motif_id_str,
                              s2_s6_open_check_cns_num_per_acr, s2_s6_target_spe2_str, s2_s6_target_spe1,
                              s2_s6_sub_open_check_CNS_overlap_ACR,
                              s2_s6_sub_open_varACRs_in_otherSpes_blastToregionInrice_shown_in_atlas,
                              s2_s7_open_check_snp_num_per_acr, input_rice_3k_snp_fl, s2_s7_target_spe1,
                              s2_s7_target_spe2_str, s2_s8_open_motif_enrich_nonsyn, s2_s8_target_spe1,
                              s2_s8_target_spe2_str,
                              s2_s8_sub_method, s2_s9_open_check_H3K27m3_under_species_specific_ACR, s2_s9_target_spe1,
                              s2_s9_target_spe2_str,
                              s2_s10_open_check_motif_enrich_on_exp_cor_ACR, s2_s10_target_spe1, s2_s10_target_ct_str,
                              s2_s10_repeat_times,
                              s2_s11_open_check_TE_on_ACR, s2_s11_target_spe1, s2_s11_target_spe2_str,
                              s2_s11_sub_open_overlap_enrich_TE, s2_s11_sub_open_motif_enrich_target_TEfam,
                              s2_s11_sub_target_CT_TEfam_str, s2_s11_sub_repeat_times,
                              s2_s11_sub_open_DNAmethylation_target_TEfam
                              )

if step03 == 'yes':
    step03_summarize_overview_H3K27me3(input_output_dir, input_all_celltype_acr_dir, input_all_spe_gene_gff_dir,
                                       input_target_organ_H3K27_peak_all_dir, input_rice_all_organ_H3K27_fl_dir,
                                       input_all_spe_PRE_motif_fimo_dir, input_all_spe_cns_gff_fl_dir,
                                       input_total_rice_atlas_acr_fl,
                                       s3_s1_open_iden_acr_to_gene_cate,
                                       s3_s2_open_summarize_overview_cross_species,
                                       s3_s2_extend_H3K27me3_flank_bp, s3_s2_target_spe1, s3_s2_target_spe2_str,
                                       s3_s3_open_check_conserved_H3K27me3_broad_ACR_in_other_tissue,
                                       s3_s3_target_organ_str,
                                       s3_s4_open_check_motif_freq_in_diff_cate_acr,
                                       s3_s5_open_check_cns_overlap_ACR, s3_s6_open_check_overlap_atlas_acr,
                                       s3_s7_open_check_target_motif_in_ACR, s3_s7_target_motifID)


if step04 == 'yes':
    step04_intersect_three_species(input_output_dir, s0_s1_target_spe1,
                                   s4_s1_open_find_candidate_plot_acr_H3K27me3,
                                   s4_s2_open_find_candidate_plot_five_spe_acr,
                                   s4_s3_open_find_candidate_plot_target_celltype, s4_s3_target_spe2_str,
                                   s4_s3_target_celltype_str)
