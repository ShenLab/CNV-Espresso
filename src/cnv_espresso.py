from __future__ import division
import os
import argparse
import glob
import function as func
import operator 
import fileinput
import pdb
import sys
import time
import itertools
from datetime import datetime
#from ParameterEstimation import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
from windows import *
from normalization import *

# Functions
def fetch_norm_rd(sampleID, sample_norm_file):
    df = pd.read_table(sample_norm_file,low_memory=False,header=None, sep='\t',
                       names=['chr', 'start','end','GC','RD',sampleID])
    return df[sampleID].to_frame()

def fetch_sampleID_from_filepath(filepath):
    filepath,tempfilename = os.path.split(filepath[0])
    sampleID = tempfilename.replace(".cov.bed.norm.gz","")
    return sampleID


# Main functions
def windows(args):
    print('Building windows.bed ...')
    output_file = args.output + '/windows.bed'
    build_windows(args.target, args.ref, output_file)

def normalization(args): 
    windows_file = str(args.windows)
    debug_flag   = args.debug
    
    input_file_list = []
    if args.input and not args.input_list:
        input_file = str(args.input)
        input_file_list.append([input_file])
    elif args.input_list:
        input_file_list = func.fileToList(args.input_list)
    else:
        print("You must input an input file or file list.")
        exit(0)

    if args.output:
        output_dir   = str(args.output)
        output_dir = output_dir 
    else:
        output_dir = os.getcwd() + '/norm'

    gc_normalization(windows_file, input_file_list, output_dir, debug_flag)   

def reference_selection(args):
    project_path   = str(args.project_path) 
    norm_list_file = str(args.norm_list)
    num_ref        = int(args.num_ref)
    corr_threshold = float(args.corr_threshold)
    
    num = 0
    combined_list = []
    sampleID_list = []
    
    # Import samples
    sample_norm_rd_list = func.fileToList_tab(norm_list_file)
    for sample_rd_path in sample_norm_rd_list:
        sampleID = fetch_sampleID_from_filepath(sample_rd_path)
        sample_rd_file = sample_rd_path[0]
        num += 1
        if num % 10 == 0:
            func.showDateTime('\t')
            print("Importing No.%d sample:%s from %s"%(num, sampleID, sample_rd_file))
        rd_df = fetch_norm_rd(sampleID, sample_rd_file)
        combined_list.append(rd_df.to_numpy())
        sampleID_list.append(sampleID)

    combined_np_array = np.hstack(combined_list) # convert to nparray to accelerate the speed.
    combined_df = pd.DataFrame(combined_np_array,columns=sampleID_list)

    # Calculate correlation coefficient for read depth matrix
    print("Calculate correlation coefficient for read depth matrix ...")
    corrMatrix = pd.DataFrame(np.corrcoef(combined_np_array, rowvar=False), columns=sampleID_list)
    sampleID_str = '\t'.join([sampleID for sampleID in sampleID_list])

    corr_matrix_file = project_path+'/correlation_matrix.txt'
    np.savetxt(corr_matrix_file,corrMatrix, delimiter="\t", header=sampleID_str, comments='')

    # Select reference samples
    for case_sampleID in sampleID_list:
        ref_sample_list = []

        print("For case sampleID:", case_sampleID)
        ref_sample_df = corrMatrix[case_sampleID].sort_values(ascending=False)
        ref_sample_size = min(num_ref,len(corrMatrix))
        for i in range(1,ref_sample_size):
            ref_sampleID = sampleID_list[ref_sample_df.index[i]]
            ref_sample_corr = ref_sample_df.iloc[i]
            if ref_sample_corr >= corr_threshold:
                ref_sample_list.append([ref_sampleID, ref_sample_corr])

        output_ref_file = project_path + '/ref_samples/'+case_sampleID+'.ref.samples.txt'
        func.output_to_file(ref_sample_list, output_ref_file)

parser = argparse.ArgumentParser(prog='cnv_espresso', description='Validate CNVs in silico')
subparsers = parser.add_subparsers()

## Transform target probe file to windows.bed and annotate GC content.
win_parser = subparsers.add_parser('windows', help="Build windows.bed file")
win_parser.add_argument('--target', required=True, help='Target probe file')
win_parser.add_argument('--ref', required=True, help='Reference file')
win_parser.add_argument('--output', required=True, help='Directory for windows.bed file')
win_parser.set_defaults(func=windows)

#Normalize read depth signal
svd_parser = subparsers.add_parser('normalization', help="GC correction, zscore by negative distribution for a given sample")
svd_parser.add_argument('--windows', required=True, help='Please input the target information including GC content')
svd_parser.add_argument('--input', required=False, help='Please input a read depth file for a given sample')
svd_parser.add_argument('--input_list', required=False, help='Please input a read depth file list for a given batch of samples')
svd_parser.add_argument('--output', required=False, help='Output folder for normalized read depth files')
svd_parser.add_argument('--debug', required=False, default=False,  help='Output folder for normalized read depth files')
svd_parser.set_defaults(func=normalization)

#Select reference samples
ref_parser = subparsers.add_parser('reference', help="Calculate the correlation matrix and select references")
ref_parser.add_argument('--project_path', required=True, help='Project folder')
ref_parser.add_argument('--norm_list', required=True, help='Normlized read depth file list')
ref_parser.add_argument('--num_ref', required=False, default=100, help='Number of reference samples')
ref_parser.add_argument('--corr_threshold', required=False, default=70, help='The mininum Pearson correlation threshold for reference samples')
ref_parser.set_defaults(func=reference_selection)

## Filter matrix by GC content, mapping ability and exon length
#svd_parser = subparsers.add_parser('filter', help="Filter matrix by GC content, mapping ability and exon length")
#svd_parser.add_argument('--rpkm_matrix', required=True, help='Matrix of RPKM values')
#svd_parser.add_argument('--ref_file', required=False, help='Reference file for the calculation of GC percentage')
#svd_parser.add_argument('--map_file', required=False, help='Mapping ability file.')
#svd_parser.add_argument('--filter_params', required=True, help='Parameters of filtering')
#svd_parser.add_argument('--output', required=False, help='Filtered matrix')
#svd_parser.set_defaults(func=filter_rpkm)
#
##SVD
#svd_parser = subparsers.add_parser('svd', help="SVD")
#svd_parser.add_argument('--rpkm_matrix', required=True, help='')
#svd_parser.add_argument('--output', required=False, help='')
## svd_parser.add_argument('--svd', type=int, required=True, help='Number of components to remove')
#svd_parser.set_defaults(func=svd)
#
##CNV discover
#cnv_parser = subparsers.add_parser('discover', help="Run HMM to discover CNVs")
#cnv_parser.add_argument('--params', required=True, help='Parameters used by HMM')
#cnv_parser.add_argument('--rpkm_matrix', required=True, help='RPKM matrix.')
#cnv_parser.add_argument('--mode',required=True, default='SVD', help='Data normalization by SVD or baseline mode.')
#cnv_parser.add_argument('--output', required=True, help='Output file.')
#cnv_parser.add_argument('--sample', required=False, default='', help='Optionally, users can choose one sample to run.')
#cnv_parser.add_argument('--vcf', required=False, default='', help='Optionally, users can input snp information by specifing a vcf file')
#cnv_parser.add_argument('--hetsnp', required=False, default=False)
##cnv_parser.add_argument('--no-hetsnp', dest='hetsnp', action='store_true')
##cnv_parser.set_defaults(hetsnp=True)
#cnv_parser.add_argument('--tagsnp', required=False, default=False)
##cnv_parser.add_argument('--no-tagsnp', dest='tagsnp', action='store_true')
##cnv_parser.set_defaults(tagsnp=True)
#cnv_parser.add_argument('--tagsnp_file',required = False, help='TagSNP file location.')
#cnv_parser.set_defaults(func=discover)
#
##Merge results
#cnv_parser = subparsers.add_parser('merge', help="Merge results from different methods")
#cnv_parser.add_argument('--datafile_svd', required=True)
#cnv_parser.add_argument('--datafile_dis', required=True)
#cnv_parser.add_argument('--output', required=True)
#cnv_parser.set_defaults(func=merge_results)
#
args = parser.parse_args()
args.func(args)
