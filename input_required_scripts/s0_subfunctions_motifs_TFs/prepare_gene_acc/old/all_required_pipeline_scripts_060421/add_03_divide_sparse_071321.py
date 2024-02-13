#!/usr/bin/env python

##this script is to divide the sparse into different parts
import re
import glob
import sys
import subprocess
import os

input_acc_sparse = sys.argv[1]
input_gba_sparse = sys.argv[2]
input_output_dir = sys.argv[3]

target_sample_string = sys.argv[4]

target_sample_list = target_sample_string.split(',')

def divide_sparse (input_acc_sparse,input_gba_sparse,target_sample_list,input_output_dir):

    ##for the acc sparse
    store_temp_cicero_acc_sparse_dir = input_output_dir + '/store_temp_cicero_acc_sparse_dir'
    if not os.path.exists(store_temp_cicero_acc_sparse_dir):
        os.makedirs(store_temp_cicero_acc_sparse_dir)

    for eachsample in target_sample_list:

        sparse_file = open (input_acc_sparse,'r')
        opt_sparse_file = open (store_temp_cicero_acc_sparse_dir + '/' + eachsample + '.sparse','w')

        for eachline in sparse_file:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if eachsample in col[1]:
                opt_sparse_file.write(eachline + '\n')

        sparse_file.close()
        opt_sparse_file.close()

    ##for the gba sparse
    store_temp_gba_sparse_dir = input_output_dir + '/store_temp_gba_sparse_dir'
    if not os.path.exists(store_temp_gba_sparse_dir):
        os.makedirs(store_temp_gba_sparse_dir)

    for eachsample in target_sample_list:

        sparse_file = open(input_gba_sparse, 'r')
        opt_sparse_file = open(store_temp_gba_sparse_dir + '/' + eachsample + '.sparse', 'w')

        for eachline in sparse_file:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if eachsample in col[1]:
                opt_sparse_file.write(eachline + '\n')

        sparse_file.close()
        opt_sparse_file.close()

divide_sparse (input_acc_sparse,input_gba_sparse,target_sample_list,input_output_dir)














