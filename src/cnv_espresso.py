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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
from windows import *
from normalization import * 
from reference import *
from images import *

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

def reference(args):
    project_path   = str(args.project_path) 
    norm_list_file = str(args.norm_list)
    num_ref        = int(args.num_ref)
    corr_threshold = float(args.corr_threshold)
    print("new")
    reference_selection(project_path, norm_list_file, num_ref, corr_threshold)

def images(args):
    RD_norm_dir     = args.rd_norm_dir
    ref_samples_dir = args.ref_dir
    cnv_file        = args.cnv_list
    output_path     = args.output
    corr_threshold  = float(args.corr_threshold)
    flanking        = args.flanking
    split_img       = args.split
    try:
        sge_task_id = int(args.specific)
    except:
        sge_task_id = 'all'

    generate_images(RD_norm_dir, ref_samples_dir, cnv_file, output_path, corr_threshold, flanking, split_img, sge_task_id)

# Main function
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
ref_parser.add_argument('--project_dir', required=True, help='Project folder')
ref_parser.add_argument('--norm_list', required=True, help='Normlized read depth file list')
ref_parser.add_argument('--num_ref', required=False, default=100, help='Max number of reference samples')
ref_parser.add_argument('--corr_threshold', required=False, default=-1, help='The mininum Pearson correlation threshold for reference samples')
ref_parser.set_defaults(func=reference)

#Generate images
img_parser = subparsers.add_parser('images', help="Encode CNV predictions into images")
img_parser.add_argument('--rd_norm_dir', required=True, help='The folder for normalized read depth files')
img_parser.add_argument('--ref_dir', required=True, help='The folder for reference samples')
img_parser.add_argument('--cnv_list', required=True, help='Please input a CNV prediction list')
img_parser.add_argument('--output', required=True, help='Output folder for images')
img_parser.add_argument('--corr_threshold', required=False, default=0.7, help='The folder for normalized read depth files')
img_parser.add_argument('--flanking', required=False, default=False, help='The folder for normalized read depth files')
img_parser.add_argument('--split', required=False, default=False, help='Generate split sliding window images for CNVs')
img_parser.add_argument('--specific', required=False, default=False, help='Generate ONE image for a specific CNV in the list file')
img_parser.set_defaults(func=images)


args = parser.parse_args()
args.func(args)
