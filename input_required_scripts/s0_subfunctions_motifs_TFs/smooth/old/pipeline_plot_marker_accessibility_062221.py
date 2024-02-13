#!/usr/bin/env python

##updating 080921 add an argument to decide to plot the known markers or not
##since we calculate the gene acc the filtered cells would be less than the total number of meta since a lots of them are removed from activity filtration
##updating 071521 add more arguments
##this script we will plot the marker acc
import re
import glob
import sys
import subprocess
import os


input_required_script_dir = sys.argv[1]

input_sparse_fl = sys.argv[2]
##/scratch/hy17471/rice_mPing_scATAC_seq_050521/07_2_calculate_gene_accessibility_060321/output_dir_12kcells_1stround_ver/03_merge_acc_dir/All.GAadjusted.sparse

input_meta_fl = sys.argv[3]
##/scratch/hy17471/rice_mPing_scATAC_seq_050521/pre_ver/06_pipe_peakcalling_peaksparse_clustering_1st_060121/output_dir_notrm_12kcells_061221/03_store_clustering_res0.25_dir/regmodel2/All.LouvainCluster.txt

input_svd_fl = sys.argv[4]
##/scratch/hy17471/rice_mPing_scATAC_seq_050521/pre_ver/06_pipe_peakcalling_peaksparse_clustering_1st_060121/output_dir_notrm_12kcells_061221/03_store_clustering_res0.25_dir/regmodel2/All.SVD.txt

input_marker_fl = sys.argv[5]
##markers_leaf.txt

input_config_fl = sys.argv[6]
##ipt_config_fl.txt

input_output_dir = sys.argv[7] ##output_dir


def Step01_de_novo_plot_markers (input_required_script_dir,input_sparse_fl,input_meta_fl,input_svd_fl,input_marker_fl,input_core_num,
                                 target_cluster,plot_each_CT,input_output_dir,onlysmooth_marker,run_denovo,run_known):

    ##create dir
    step01_firstround_dir = input_output_dir + '/step01_firstround_create_obj_dir'
    if not os.path.exists(step01_firstround_dir):
        os.makedirs(step01_firstround_dir)

    ##the first round will run all the markers at the first time
    ##if we obtain file from the first round we will consider use the opt to generate markers as we want
    plot_marker_script = input_required_script_dir + '/Step01_plot_marker_accessibility_062221.R'
    function_script = input_required_script_dir + '/Step01_functions.plot_marker_accessibility.R'

    cmd = 'Rscript ' + plot_marker_script + \
          ' ' + input_meta_fl + \
          ' ' + input_sparse_fl + \
          ' ' + input_svd_fl + \
          ' ' + input_marker_fl + \
          ' ' + input_core_num + \
          ' ' + target_cluster + \
          ' ' + plot_each_CT + \
          ' ' + step01_firstround_dir + \
          ' ' + function_script + \
          ' ' + onlysmooth_marker + \
          ' ' + run_denovo + \
          ' ' + run_known
    subprocess.call(cmd,shell=True)


def Step02_only_plot_provided_markers (input_required_script_dir,plot_each_CT,input_marker_fl,input_output_dir,lim):


    ##create dir
    step02_onlyplot_provided_marker_dir = input_output_dir + '/step02_onlyplot_provided_marker_dir'
    if not os.path.exists(step02_onlyplot_provided_marker_dir):
        os.makedirs(step02_onlyplot_provided_marker_dir)

    ##updating 071521
    ##find the activity.rds
    all_files = glob.glob(input_output_dir + '/step01_firstround_create_obj_dir/*impute.activity.rds')
    all_gene_impute_acc = all_files[0]

    #all_gene_impute_acc = input_output_dir + '/step01_firstround_create_obj_dir/opt_allgenes_impute.activity.rds'
    GAobj_rds = input_output_dir + '/step01_firstround_create_obj_dir/GAobj.rds'

    step02_script = input_required_script_dir + '/Step02_only_plot_provided_marker_accessibility_062221.R'

    cmd = 'Rscript ' + step02_script + \
          ' ' + all_gene_impute_acc + \
          ' ' + plot_each_CT + \
          ' ' + GAobj_rds + \
          ' ' + input_marker_fl + \
          ' ' + step02_onlyplot_provided_marker_dir + \
          ' ' + lim
    subprocess.call(cmd,shell=True)



##open the configure file to check the parameter settings
store_target_parameter_dic = {}
with open (input_config_fl,'r') as ipt:
    for eachline in ipt:
        eachline = eachline.strip('\n')
        if not eachline.startswith('#'):
            col = eachline.strip().split('=')
            store_target_parameter_dic[col[0]] = col[1]

input_core_num = store_target_parameter_dic['core']
target_cluster = store_target_parameter_dic['target_cluster']
plot_each_CT = store_target_parameter_dic['plot_each_CT']
step01 = store_target_parameter_dic['step01'] ##step01 is yes or no
step02 = store_target_parameter_dic['step02'] ##step02 is yes or no
lim = store_target_parameter_dic['lim']
onlysmooth_marker = store_target_parameter_dic['onlysmooth_marker']
run_denovo = store_target_parameter_dic['run_denovo']
run_known = store_target_parameter_dic['run_known']


if step01 == 'yes':
    Step01_de_novo_plot_markers (input_required_script_dir,input_sparse_fl,input_meta_fl,input_svd_fl,input_marker_fl,input_core_num,
                                     target_cluster,plot_each_CT,input_output_dir,onlysmooth_marker,run_denovo,run_known)

if step02 == 'yes':
    Step02_only_plot_provided_markers(input_required_script_dir, plot_each_CT, input_marker_fl, input_output_dir, lim)





