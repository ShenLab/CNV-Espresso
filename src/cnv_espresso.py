from __future__ import division
import os
import argparse
import glob
import function as df
#from DataManager import *
#from hmm.Model import *
#from hmm.ModelParams import *
import operator 
import numpy as np
#from VCFReader import *
from ParameterEstimation import *
import fileinput
import pdb
import sys
import time

def normalization(args): 
    windows_file = str(args.windows)
    input_file_list = []
    if args.input and not args.input_list:
        input_file = str(args.input)
        input_file_list.append(input_file)

    elif args.input_list:
        input_file_list = df.fileToList(args.input_list)
    else:
        print("You must input an input file or file list.")
        exit(0)

    if args.output:
        output_dir   = str(args.output)
        output_dir = output_dir + '/'
    else:
        output_dir = os.getcwd() 

    print('Loading windows.bed ...')
    windows_dict = df.loadWindows(windows_file)
    windows_chr  = windows_dict['chr']
    windows_start= windows_dict['start']
    windows_stop = windows_dict['stop']
    windwos_gc   = windows_dict['gc']
    windows_mappability = windows_dict['mappability']

    for input_file in input_file_list:
        input_file = input_file
        (input_dir,input_name) = os.path.split(input_file)
        (input_filename,input_extension) = os.path.splitext(input_name)

        output_file = output_dir + input_filename + '.norm'
        output_parameter_file = output_dir + input_filename + '.nb.parm'

        print('Loading %s ...'% input_file)
        input_sample_dict = df.loadRD(input_file)
        input_sample_chr   = input_sample_dict['chr']
        input_sample_start = input_sample_dict['start']
        input_sample_stop  = input_sample_dict['stop']
        input_sample_rd    = input_sample_dict['RD']

        if not np.array_equal(input_sample_chr, windows_chr) or \
                not np.array_equal(input_sample_start, windows_start) or \
                not np.array_equal(input_sample_stop, windows_stop):
                    print("[Error] The windows file is not consisitant with input sample RD file.")
                    sys.exit(0)

        GC_percentage = np.round(windwos_gc*100)
        GC_index = {}
        for ind in range(len(GC_percentage)):
            gc = int(GC_percentage[ind])

            if gc in GC_index:
                GC_index[gc].append(ind)
            else:
                GC_index[gc] = [ind]

        print( 'Normalizing by GC percentage...')
        corrected_rd = np.zeros(len(input_sample_rd), dtype=np.float)
        overall_median = np.median(GC_percentage)
        for gc in GC_index.keys():
            t_ind = GC_index[gc]
            t_median = np.median(input_sample_rd[t_ind])
            if t_median == 0:
                print( 'WARNING. Median == 0, GC: %d' %(gc))
                corrected_rd[t_ind] = 0
            else:
                corrected_rd[t_ind] = np.round(input_sample_rd[t_ind] * overall_median / t_median)


        # Estimate NB distribution  parameters
        parameterLoader = ParameterEstimation(corrected_rd)
        parameterList = parameterLoader.fit(corrected_rd,0.01,0.99)#TODO remove sex chromsomes and estimate them seperately
        print("Estimated Paramters: ", parameterList)
        mu = parameterList[0]
        fi = parameterList[1]
        sigma = (mu + (mu**2)*fi)**0.5
        corrected_zscore = (corrected_rd - mu)/sigma

        # Output
        output_ndarray = np.transpose(np.array([input_sample_chr, input_sample_start, input_sample_stop, \
                GC_percentage, input_sample_rd, corrected_rd, corrected_zscore]))
        print( 'Saving normalized read depth to file %s'%output_file)
        df.output_to_file(output_ndarray, output_file)
        df.output_to_file(parameterList, output_parameter_file)


parser = argparse.ArgumentParser(prog='cnv_espresso', description='Validate CNVs in silico')
subparsers = parser.add_subparsers()

##BAM List -> RPKM
#svd_parser = subparsers.add_parser('rpkm', help="Create RPKM matrix from a BAM list")
#svd_parser.add_argument('--target', required=True, help='Target definition file')
#svd_parser.add_argument('--maq', required=False, type=int, default=20, help='MAQ threshold')
#svd_parser.add_argument('--input', required=True, help='BAM file list, each line for each sample')
#svd_parser.add_argument('--output', required=True, help='Directory for RPKM files')
#svd_parser.set_defaults(func=bamlist2RPKM)
#
#normalize read depth signal
svd_parser = subparsers.add_parser('normalization', help="GC correction, zscore by negative distribution for a given sample")
svd_parser.add_argument('--windows', required=True, help='Please input the target information including GC content')
svd_parser.add_argument('--input', required=False, help='Please input a read depth file for a given sample')
svd_parser.add_argument('--input_list', required=False, help='Please input a read depth file list for a given batch of samples')
svd_parser.add_argument('--output', required=False, help='Output folder for normalized read depth files')
svd_parser.set_defaults(func=normalization)

##RPKM files -> Matrix
#svd_parser = subparsers.add_parser('merge_rpkm', help="Merge RPKM files to a matrix")
#svd_parser.add_argument('--rpkm_dir', required=True, help='RPKM files')
#svd_parser.add_argument('--target', required=True, help='Target definition file')
#svd_parser.add_argument('--output', required=False, help='Matrix file')
#svd_parser.set_defaults(func=RPKM2Matrix)
#
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
