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
from windows import *
from normalization import * 
from reference import *
from images import *
from images_human_view import *
from train import *
from predict import *

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
    project_dir    = str(args.project_dir) 
    norm_list_file = str(args.norm_list)
    num_ref        = int(args.num_ref)
    corr_threshold = float(args.corr_threshold)
    reference_selection(project_dir, norm_list_file, num_ref, corr_threshold)

def images(args):
    RD_norm_dir     = args.rd_norm_dir
    ref_samples_dir = args.ref_dir
    cnv_file        = args.cnv_list
    output_path     = args.output
    corr_threshold  = float(args.corr_threshold)
    flanking        = args.flanking
    split_img       = args.split
    job_start       = args.start
    overwrite_img   = args.overwrite_img
    try:
        sge_task_id = int(args.specific)
    except:
        sge_task_id = 'all'

    generate_images(RD_norm_dir, ref_samples_dir, cnv_file, output_path, corr_threshold, \
                    flanking, split_img, sge_task_id, job_start, overwrite_img)

def images_human_view(args):
    RD_norm_dir     = args.rd_norm_dir
    ref_samples_dir = args.ref_dir
    cnv_file        = args.cnv_list
    vcf_file        = args.vcf_file
    output_path     = args.output
    suffix          = args.suffix
    corr_threshold  = float(args.corr_threshold)
    flanking        = args.flanking
    split_img       = args.split
    trio            = args.trio
    pedigree_file   = args.pedigree
    overwrite_img   = args.overwrite_img
    try:
        sge_task_id = int(args.specific)
    except:
        sge_task_id = 'all'

    generate_images_human_view(RD_norm_dir, ref_samples_dir, cnv_file, vcf_file, output_path, suffix, corr_threshold, flanking, \
                               split_img, sge_task_id, trio, pedigree_file, overwrite_img)

def train(args):
    true_del_file    = args.true_del
    true_dup_file    = args.true_dup
    false_del_file   = args.false_del
    false_dup_file   = args.false_dup
    use_gpu          = args.use_gpu
    batch_size       = args.batch_size
    epochs           = args.epochs
    output_model_dir = args.output

    cnn_train(true_del_file, true_dup_file, false_del_file, false_dup_file, use_gpu, batch_size, epochs, output_model_dir)

def predict(args):
    cnv_file    = args.cnv_list
    model_file  = args.model
    output_file = args.output
    use_gpu     = args.use_gpu
    cnn_prediction(cnv_file, model_file, use_gpu, output_file)

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
img_parser.add_argument('--start', required=False, default=False, help='The number from which image is generated')
img_parser.add_argument('--overwrite_img', required=False, default=True, help='Overwrite the current image if it exists.')
img_parser.set_defaults(func=images)

#Generate images for human view
img_parser = subparsers.add_parser('images_human_view', help="Encode CNV predictions into images for human review")
img_parser.add_argument('--rd_norm_dir', required=True, help='The folder for normalized read depth files')
img_parser.add_argument('--ref_dir', required=True, help='The folder for reference samples')
img_parser.add_argument('--cnv_list', required=True, help='Please input a CNV prediction list')
img_parser.add_argument('--vcf_file', required=False, help='Please input a VCF file which include SNV information for this region')
img_parser.add_argument('--output', required=True, help='Output folder for images')
img_parser.add_argument('--suffix', required=False, default=None, help='add the suffix to the end of output file name')
img_parser.add_argument('--corr_threshold', required=False, default=0.7, help='The folder for normalized read depth files')
img_parser.add_argument('--flanking', required=False, default=False, help='The folder for normalized read depth files')
img_parser.add_argument('--split', required=False, default=False, help='Generate split sliding window images for CNVs')
img_parser.add_argument('--specific', required=False, default=False, help='Generate ONE image for a specific CNV in the list file')
img_parser.add_argument('--trio', required=False, default=False, help='Trio mode or single sample mode')
img_parser.add_argument('--pedigree', required=False, default=False, help='If in trio mode, a pedigree file is required')
img_parser.add_argument('--overwrite_img', required=False, default=True, help='Overwrite the current image if it exists.')
img_parser.set_defaults(func=images_human_view)

#Train the CNN model
cnn_parser = subparsers.add_parser('train', help="Train the CNN model from scratch")
cnn_parser.add_argument('--true_del', required=True, help='Please input a bunch of true deletions with image path for training the model')
cnn_parser.add_argument('--true_dup', required=True, help='Please input a bunch of true duplications with image path for training the model')
cnn_parser.add_argument('--false_del', required=True, help='Please input a bunch of false deletions with image path for training the model')
cnn_parser.add_argument('--false_dup', required=True, help='Please input a bunch of false duplications with image path for training the model')
cnn_parser.add_argument('--use_gpu', required=False, default=True, help='Using GPU or CPU to train the model')
cnn_parser.add_argument('--batch_size', required=False, default=32, help='Please input the batch size for CNN training')
cnn_parser.add_argument('--epochs', required=False, default=20, help='Please input epochs for CNN training')
cnn_parser.add_argument('--output', required=True, help='Output directory for a trained CNN model')
cnn_parser.set_defaults(func=train)

#Prediction
pred_parser = subparsers.add_parser('predict', help="Validate CNV predictions by CNN")
pred_parser.add_argument('--cnv_list', required=True, help='Please input the CNV list for validation')
pred_parser.add_argument('--model', required=True, help='Please input a trained CNN model file')
pred_parser.add_argument('--use_gpu', required=False, default=False, help='Predict CNVs by using GPU or CPU')
pred_parser.add_argument('--output', required=True, help='Output CNV validation results')
pred_parser.set_defaults(func=predict)

args = parser.parse_args()
args.func(args)
