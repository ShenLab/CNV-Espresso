#!/usr/bin/env python
# coding: utf-8
# # Generate images for both human eyes and in silico
import function as func
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pdb
import os
import bz2
import math
import pysam
import sys
import glob
import copy

## Variables
SAMPLE      = func.global_variables()['SAMPLE']
CNV_INTERVAL= func.global_variables()['CNV_INTERVAL']
CNV_CHR     = func.global_variables()['CNV_CHR']
CNV_START   = func.global_variables()['CNV_START']
CNV_END     = func.global_variables()['CNV_END']
CNV_TYPE    = func.global_variables()['CNV_TYPE']
NUM_TARGETS = func.global_variables()['NUM_TARGETS']
CNV_LABEL   = func.global_variables()['CNV_LABEL']

target_group = 3 # the number of targets per group
color_del = (0,1,0) #green
color_dup = (1,0,0) #red

## Input
cnv_file        = sys.argv[1] 
output_path     = sys.argv[2]
#output_prefix   = sys.argv[3]

## Output files and folders
output_EntireCNV_image_dir  = output_path + '/entire_cnvs/'
output_SplitCNV_image_dir   = output_path + '/split_cnvs/'
os.makedirs(output_EntireCNV_image_dir,  exist_ok=True)
os.makedirs(output_SplitCNV_image_dir,   exist_ok=True)

## Functions
def fetch_colName(keyWord_list, colName_list):
    keyWord = [i for i in colName_list if i in keyWord_list]
    if len(keyWord) == 1:
        return keyWord[0]
    else:
        print("[Error] Please check the keyWord: "%';'.join([s for s in keyWord]))
        #pdb.set_trace()
        return None

## CNV info
print("Import cnv data ...")
cnv_data_df = pd.read_table(cnv_file, header=0)

## Processing
## Parse header
print("Parsing header ...")
cnv_data_header = cnv_data_df.columns.tolist()
cnv_data_header = [col_name.upper() for col_name in cnv_data_header] # upper all colnames

col_sampleID  = cnv_data_header.index(fetch_colName(cnv_data_header,SAMPLE))
try:
    col_cnv_interval = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_INTERVAL))
except:
    col_cnv_chr   = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_CHR))
    col_cnv_start = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_START))
    col_cnv_end   = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_END))
col_cnv_type  = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_TYPE))

try:
    col_cnv_num_targets = cnv_data_header.index(fetch_colName(cnv_data_header,NUM_TARGETS))
    col_cnv_label = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_LABEL))
    col_cnv_canoes= cnv_data_header.index(fetch_colName(cnv_data_header,['CANOES','CANOES_RT']))
    col_cnv_xhmm  = cnv_data_header.index(fetch_colName(cnv_data_header,['XHMM','XHMM_RT']))
    col_cnv_clamms= cnv_data_header.index(fetch_colName(cnv_data_header,['CLAMMS','CLAMMS_RT']))
    col_cnv_num_carriers = cnv_data_header.index(fetch_colName(cnv_data_header,['Num_Carriers(inGivenCohort)']))
    col_cnv_num_wins  = cnv_data_header.index(fetch_colName(cnv_data_header,['Num_Targets_Wins','Num_Targets','Num_Wins']))
except:
    pass
## Checking results

### add image path into the CNV_info file
output_df = copy.deepcopy(cnv_data_df)
#output_df.loc[:, 'image_path'] = ["-"]*len(cnv_data_df) 
output_df.loc[:,'entire_cnv_path'] = ['-']*len(output_df)
output_df.loc[:,'split_cnv_path']  = ['-']*len(output_df)

print("Checking results ...")
miss_num = 0
for index, row in cnv_data_df.iterrows(): 
    row      = cnv_data_df.iloc[index]
    sampleID = row[col_sampleID]
    try:
        cnv_interval = row[col_cnv_interval]
        cnv_chr,cnv_start,cnv_end = func.parseInterval(cnv_interval)
    except:
        cnv_chr   = row[col_cnv_chr]
        cnv_start = int(row[col_cnv_start])
        cnv_end   = int(row[col_cnv_end])
    cnv_type = row[col_cnv_type]

    if cnv_type == 1:
        cnv_type = "DEL"
    elif cnv_type == 3:
        cnv_type = "DUP"
    else:
        pass

    try:
        cnv_num_targets   = row[col_cnv_num_targets]
        cnv_label         = row[col_cnv_label]
        cnv_canoes        = str(row[col_cnv_canoes])
        cnv_xhmm          = str(row[col_cnv_xhmm])
        cnv_clamms        = str(row[col_cnv_clamms])
        case_sample_color = color_del if cnv_type == 'DEL' else color_dup
        cnv_num_carriers  = str(row[col_cnv_num_carriers])
        cnv_num_wins      = str(row[col_cnv_num_wins])
    except:
        cnv_num_targets   = 'NA' 
        cnv_label         = 'NA'
        cnv_canoes        = 'NA'
        cnv_xhmm          = 'NA'
        cnv_clamms        = 'NA'
        cnv_num_carriers  = 'NA'
        cnv_num_wins      = 'NA'
    cnv_label_str = "True" if cnv_label == 1 else "False"
    #whole cnv image
    image_file = str(index+1)+"_"+sampleID+"_"+str(cnv_chr)+"_"+str(cnv_start)+"_"+str(cnv_end) + \
                 "_"+str(cnv_num_targets)+"tgs_"+str(cnv_num_wins)+"wins_"+ \
                 cnv_label_str+"_"+cnv_type+".png"

    ## Core. check whether the image file exists and annotate the path to CNV file
    ### entire CNVs
    if os.path.exists(output_EntireCNV_image_dir+image_file):
        full_path = output_EntireCNV_image_dir+image_file
    else:
        image_prefix = sampleID+"_"+str(cnv_chr)+"_"+ \
                        str(cnv_start)+"_"+str(cnv_end)+"_*"+cnv_type
        image_likely_file = glob.glob(output_EntireCNV_image_dir+'*_'+image_prefix+'*.png')
        if len(image_likely_file) == 1:
            #print("Find related file: %s "%(image_likely_file))
            full_path = image_likely_file[0]
        else:
            miss_num += 1
            print("%d [%d|%d] Missing results: %s "%(miss_num, index+1, len(cnv_data_df), image_prefix))
            pdb.set_trace()    
    full_path = full_path.replace("//","/")        
    output_df.loc[index, 'entire_cnv_path'] = full_path 
    
    ### split CNVs
    if os.path.exists(output_SplitCNV_image_dir+image_file):
        split_path = output_SplitCNV_image_dir+image_file
    else:
        image_prefix = sampleID+"_"+str(cnv_chr)+"_"+ \
                        str(cnv_start)+"_"+str(cnv_end)+"_*"+cnv_type
        image_likely_file_list = glob.glob(output_SplitCNV_image_dir+'*_'+image_prefix+'*.png')
        split_cnv_path = '\n'.join([cnv_path for cnv_path in image_likely_file_list])
        split_cnv_path = split_cnv_path.replace("//","/")        
    output_df.loc[index, 'split_cnv_path'] = split_cnv_path 
    
    ## check the status
    if index % 100 == 0:
        print(index, output_df.loc[index,:])

## output results
#del_result_df = output_df[output_df["TYPE"]=='DEL']
#dup_result_df = output_df[output_df["TYPE"]=='DUP']
#del_result_df.to_csv(output_path+'/'+output_prefix+'_del_image_info.list',sep='\t',index=False)
#dup_result_df.to_csv(output_path+'/'+output_prefix+'_dup_image_info.list',sep='\t',index=False)
path, filename, file_extension = func.extractFilePathNameExtension(cnv_file)
output_file = output_path+'/'+filename+'_withImagePath.csv'
output_df.to_csv(output_file, index=False)
print("Output to file:",output_file)
