#!/bin/bash
#SBATCH --partition=schmitz_p
#SBATCH --job-name=Eseedling
#SBATCH --ntasks=1
#SBATCH --time=100:00:00
#SBATCH --mem=124G
#SBATCH --cpus-per-task=1

# set env
#cd $PBS_O_WORKDIR
#source ~/.zshrc
module load Python/3.7.4-GCCcore-8.3.0
#module load Trim_Galore/0.6.5-GCCcore-8.2.0-Java-11
#module load Perl/5.30.0-GCCcore-8.3.0
#module load FastQC/0.11.8-Java-11
#module load STAR/2.7.1a-foss-2016b
#module load SAMtools/1.10-GCC-8.3.0
#module load picard/2.16.0-Java-1.8.0_144
#module load BEDTools/2.29.2-GCC-8.3.0
#ml R/4.0.0-foss-2019b
ml R/4.2.1-foss-2020b

#cd /scratch/hy17471/soybean_scATAC_100120/pipeline_analysis_110220/add_06_geneAccessibility_UMAP_analysis_regmodel2_res01_shift_011421
#cd /scratch/hy17471/soybean_scATAC_100120/pipeline_analysis_110220/add_06_geneAccessibility_UMAP_analysis_regmodel2_res01_011321

# run
python pipeline_plot_marker_accessibility_062221.py \
input_required_script_dir \
/scratch/hy17471/rice_altlas_scATAC_seq_042021/03_2_calculate_gene_accessibility_070921/output_dir_Eseedling_ver3_Socrates_fixissue_geneVer_041722/02_call_accessiblity_normGBA_dir/genebody_accessibility_dir/Eseedling.GBaccessibility.sparse \
/scratch/hy17471/rice_altlas_scATAC_seq_042021/resources/opt_final_comb_annot_ver7.5.txt \
/scratch/hy17471/rice_altlas_scATAC_seq_042021/06_1_seperate_clustering_112421/Eseedling/output_dir_win50k_PCs30_NMF/opt_tfidf_NMF_30PCs_win50000_res1_rHM_reduced_dimensions.txt \
markers_allorgans_TFsRMDUPoriFTargetaddCommGene070623.txt \
ipt_config_fl_Eseedling_CTlevel.txt \
output_dir_win50k_PCs30_NMF_geneVer_CTlevel_v7.5
