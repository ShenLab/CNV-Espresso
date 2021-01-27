#!/usr/bin/env python
# coding: utf-8
# # Generate images for both human eyes and in silico
import function as df
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

## Variables
SAMPLE      = ['SAMPLE','sample_ID']
CNV_CHR     = ['chr', 'CHR', 'CHROMOSOME', 'chromosome']
CNV_START   = ['cnv_start', 'start', 'PRED_START', 'START']
CNV_END     = ['cnv_stop', 'stop', 'PRED_END', 'END']
CNV_TYPE    = ['cnv_type','type','TYPE','CNV', 'CNV_TYPE']
NUM_TARGETS = ['NUM_TARGETS','targets']
CNV_LABEL   = ['LABEL_VAL','label','LABEL']

target_group = 3 # the number of targets per group
color_del = (0,1,0) #green
color_dup = (1,0,0) #red

## Input
RD_norm_dir     = sys.argv[1]
ref_samples_dir = sys.argv[2]
cnv_file        = sys.argv[3] 
output_path     = sys.argv[4]
sge_task_id     = int(sys.argv[5])

## Output files and folders
output_false_del_image_dir  = output_path + '/false_del/'
output_false_del_image_splits_dir = output_path + '/false_del_splits/'

output_false_dup_image_dir  = output_path + '/false_dup/'
output_false_dup_image_splits_dir = output_path + '/false_dup_splits/'

output_true_del_image_dir  = output_path + '/true_del/'
output_true_del_image_splits_dir = output_path + '/true_del_splits/'

output_true_dup_image_dir  = output_path + '/true_dup/'
output_true_dup_image_splits_dir = output_path + '/true_dup_splits/'

os.makedirs(output_false_del_image_dir, exist_ok=True)
os.makedirs(output_false_dup_image_dir, exist_ok=True)
os.makedirs(output_true_del_image_dir,  exist_ok=True)
os.makedirs(output_true_dup_image_dir,  exist_ok=True)
os.makedirs(output_false_del_image_splits_dir, exist_ok=True)
os.makedirs(output_false_dup_image_splits_dir, exist_ok=True)
os.makedirs(output_true_del_image_splits_dir , exist_ok=True)
os.makedirs(output_true_dup_image_splits_dir , exist_ok=True)

## CNV info
cnv_data_df = pd.read_table(cnv_file, header=0)
## Functions
def fetchRDdata_byTabix(RD_norm_dir, sampleID, cnv_chr, cnv_start, cnv_end, target_group):
    # tabix RD file to fetch 
    RD_norm_file = RD_norm_dir+sampleID+'.cov.bed.norm.gz'
    if not os.path.exists(RD_norm_dir):
        print('No tabular file: %s'%RD_norm_file)
        return ["No tabular file"]
    if not os.path.exists(RD_norm_file+'.tbi'):
        pysam.tabix_index(RD_norm_file, seq_col=0, start_col=1, end_col=2) # Need to add '-p bed'
    # fetch
    f = pysam.TabixFile(RD_norm_file)
    RD_fetched_data = f.fetch(cnv_chr, int(cnv_start), int(cnv_end), parser=pysam.asTuple())
    RD_fetched_df = pd.DataFrame(data=RD_fetched_data, columns=['chr', 'start', 'end', 'GC', 'RD_raw', 'RD_norm'])
    # add a new column as target groups
    RD_fetched_df_tmp = RD_fetched_df.copy()
    RD_fetched_df_tmp.loc[:, 'target_group'] = [val for val in np.arange(1,math.ceil((len(RD_fetched_df_tmp)/target_group))+1) for i in range(target_group)][0:len(RD_fetched_df_tmp)]
    RD_fetched_df = RD_fetched_df_tmp
    del RD_fetched_df_tmp
    # change the type of columns
    RD_fetched_df[["start"]] = RD_fetched_df[["start"]].astype(int)
    RD_fetched_df[["end"]] = RD_fetched_df[["end"]].astype(int)
    RD_fetched_df[["GC"]] = RD_fetched_df[["GC"]].astype(float)
    RD_fetched_df[["RD_raw"]] = RD_fetched_df[["RD_raw"]].astype(float)
    RD_fetched_df[["RD_norm"]] = RD_fetched_df[["RD_norm"]].astype(float)
    return RD_fetched_df

def loadRefSamplesID(ref_samples_file):
    ref_samplesID_df = pd.read_table(ref_samples_file,low_memory=False,header=None, sep=' ',                              names=['sampleID', 'r2'])
    return ref_samplesID_df

def fetchRefRDdata_byTabix(ref_samples_file, cnv_chr, cnv_start, cnv_end, target_group):
    # load reference sample ID
    ref_samplesID_df = loadRefSamplesID(ref_samples_file)
    # load RD normalized data and fetch RD given the cnv region for each reference sample
    reference_RD_df = pd.DataFrame(columns=['chr', 'start', 'end', 'GC', 'RD_raw', 'RD_norm', 'sample'])
    for index, row in ref_samplesID_df.iterrows():  
        ref_sampleID = row[0]
        try:
            ref_RD_cnv_region = fetchRDdata_byTabix(RD_norm_dir, ref_sampleID, cnv_chr, cnv_start, cnv_end, target_group)
            # add a new column as sampleID
            RD_cnv_region_tmp = ref_RD_cnv_region.copy()
            RD_cnv_region_tmp.loc[:, 'sample'] = [ref_sampleID]*len(ref_RD_cnv_region)
            ref_RD_cnv_region = RD_cnv_region_tmp
            del RD_cnv_region_tmp
            # combine results
            reference_RD_df = reference_RD_df.append(ref_RD_cnv_region)
        except:
            print("    -[Error]: error in normalized reference RD file of %s in %s"%(ref_sampleID, RD_norm_dir))
            #TODO: write to log file
    return reference_RD_df

def fetch_colName(keyWord_list, colName_list):
    for keyWord in keyWord_list:
        if keyWord in colName_list:
            return keyWord
        else:
            keyWord = None
    return keyWord


# ## Processing
## Parse header
cnv_data_header = cnv_data_df.columns.tolist()
col_sampleID  = cnv_data_header.index(fetch_colName(cnv_data_header,SAMPLE))
col_cnv_chr   = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_CHR))
col_cnv_start = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_START))
col_cnv_end   = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_END))
col_cnv_type  = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_TYPE))
col_cnv_num_targets = cnv_data_header.index(fetch_colName(cnv_data_header,NUM_TARGETS))
col_cnv_label = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_LABEL))
col_cnv_canoes= cnv_data_header.index(fetch_colName(cnv_data_header,['CANOES','CANOES_RT']))
col_cnv_xhmm  = cnv_data_header.index(fetch_colName(cnv_data_header,['XHMM','XHMM_RT']))
col_cnv_clamms= cnv_data_header.index(fetch_colName(cnv_data_header,['CLAMMS','CLAMMS_RT']))
col_cnv_numCarriers = cnv_data_header.index(fetch_colName(cnv_data_header,['Num_Carriers(inGivenCohort)']))

#for index, row in cnv_data_df.iterrows(): 
#     row = next(cnv_data_df.iterrows())[1]
#     index = 1
#    if index + 1 < start_point:
#        continue
index = sge_task_id-1
row   = cnv_data_df.iloc[index]

sampleID  = row[col_sampleID]
cnv_chr   = row[col_cnv_chr]
cnv_start = np.int(row[col_cnv_start])
cnv_end   = np.int(row[col_cnv_end])
cnv_type  = row[col_cnv_type]
cnv_num_targets = row[col_cnv_num_targets]
cnv_label  = row[col_cnv_label]
cnv_canoes = str(row[col_cnv_canoes])
cnv_xhmm   = str(row[col_cnv_xhmm])
cnv_clamms = str(row[col_cnv_clamms])
case_sample_color = color_del if cnv_type == 'DEL' else color_dup
cnv_num_carriers  = str(row[col_cnv_numCarriers])

if cnv_type == 'DEL' and cnv_label == 0:
    output_image_dir = output_false_del_image_dir 
    output_image_splits_dir = output_false_del_image_splits_dir
elif cnv_type == 'DEL' and cnv_label == 1:
    output_image_dir = output_true_del_image_dir
    output_image_splits_dir = output_true_del_image_splits_dir
elif cnv_type == 'DUP' and cnv_label == 0:
    output_image_dir = output_false_dup_image_dir
    output_image_splits_dir = output_false_dup_image_splits_dir
elif cnv_type == 'DUP' and cnv_label == 1: 
    output_image_dir = output_true_dup_image_dir
    output_image_splits_dir = output_true_dup_image_splits_dir
else:
    print("cnv_label error?", cnv_label)
    pdb.set_trace()
    
cnv_label_str = "True" if cnv_label == 1 else "False"
print("[%d|%d] Illustrating: %s %s:%d-%d %s #targets:%d Label:%s"% \
       (len(cnv_data_df), index+1, sampleID, cnv_chr, cnv_start, cnv_end, cnv_type, cnv_num_targets, cnv_label_str))

## Import RD data info
print("  --Step1. Fetching RD data for case sample ...")
sample_rd_file = RD_norm_dir+sampleID+'.cov.bed.norm.gz'
if not os.path.exists(sample_rd_file):
    print("    -[Error]: error in normalized RD file of %s "%(sample_rd_file))
    #TODO: write to log file
    exit(0)
RD_cnv_region = fetchRDdata_byTabix(RD_norm_dir, sampleID, cnv_chr, cnv_start, cnv_end, target_group)

## Fetch Read depth data for reference samples in terms of CNV boundary       
print("  --Step2. Fetching RD data for reference samples ...")
ref_samples_file = ref_samples_dir+sampleID+'.ref.samples.txt.bz2'
if not os.path.exists(ref_samples_file):
    print("    -[Error]: error in reference samples related file for %s in %s"%(sampleID, ref_samples_dir))
    #TODO: write to log file
    exit(0)
reference_RD_df = fetchRefRDdata_byTabix(ref_samples_file, cnv_chr, cnv_start, cnv_end, target_group)
    
## plot whole cnv
print("  --Step3. Illustrating an image for the whole CNV ...")
title_info = sampleID+" "+str(cnv_chr)+":"+str(cnv_start)+"-"+str(cnv_end)+" "+cnv_type + \
             " "+ str((cnv_end-cnv_start)/1000) + 'kb'+ " #targets:"+str(cnv_num_targets) + \
             " #wins:" + str(len(RD_cnv_region)) + "\nCANOES:"+cnv_canoes + " XHMM:"+ \
             cnv_xhmm + " CLAMMS:"+cnv_clamms + " #Carriers:"+cnv_num_carriers

image_file = str(index+1)+"_"+sampleID+"_"+str(cnv_chr)+"_"+str(cnv_start)+"_"+str(cnv_end)+ \
             "_"+cnv_type+ "_"+str(cnv_num_targets)+"tgs_"+str(len(RD_cnv_region)) +"wins.png"

fig = plt.figure(dpi=150,figsize=(10, 7)) 
ax_rd = fig.subplots(nrows=1, ncols=1)

### plot reference samples
for sample_reader in reference_RD_df["sample"].unique():
            ref_sample_df = reference_RD_df[reference_RD_df["sample"]==sample_reader]
            ax_rd.plot((ref_sample_df["start"]+ref_sample_df["end"])/2, ref_sample_df["RD_norm"], color='grey', marker='.', linewidth=0.2)

### plot case sample
ax_rd.plot((RD_cnv_region["start"]+RD_cnv_region["end"])/2, RD_cnv_region["RD_norm"],color=case_sample_color , marker='o', linewidth=2)
ax_rd.set_title(title_info)
plt.savefig(output_image_dir+image_file)
plt.close() 

## plot split CNV for each three targets
print("  --Step4. Illustrating images for the CNV splited by each %d windows ..."%target_group)
for group_id in np.unique(RD_cnv_region['target_group']):
    ## if targets equal to required number (3 by default)
    if len(RD_cnv_region[RD_cnv_region['target_group']==group_id]) == target_group:
        title_split_info = title_info +" Group:"+ str(int(len(RD_cnv_region)/target_group))+"-"+ str(group_id)
        image_split_file = str(index+1)+"_"+sampleID+"_"+str(cnv_chr)+"_"+str(cnv_start)+"_"+str(cnv_end)+"_"+cnv_type+ "_"+str(cnv_num_targets)+                                 "tgs_"+str(len(RD_cnv_region)) +"wins_splits"+str(int(len(RD_cnv_region)/target_group))+"_"+ str(group_id) +".png"
        fig = plt.figure(dpi=150,figsize=(7, 7)) 
        ax_rd = fig.subplots(nrows=1, ncols=1)

        RD_cnv_region_split = RD_cnv_region[RD_cnv_region["target_group"]==group_id]
        reference_RD_df_split = reference_RD_df[reference_RD_df["target_group"]==group_id]
        # plot reference samples
        for sample_reader in reference_RD_df_split["sample"].unique():
            ref_sample_df = reference_RD_df_split[reference_RD_df_split["sample"]==sample_reader]
            ax_rd.plot((ref_sample_df["start"]+ref_sample_df["end"])/2, ref_sample_df["RD_norm"], color = 'grey', marker='.', linewidth=0.2)
        # plot case sample
        ax_rd.plot((RD_cnv_region_split["start"]+RD_cnv_region_split["end"])/2, RD_cnv_region_split["RD_norm"], color=case_sample_color, marker='o', linewidth=2)
        ax_rd.set_title(title_split_info)
        plt.savefig(output_image_splits_dir+image_split_file)
        plt.close()
print("  --[Done]. Images have output to %s and %s."%(output_image_dir+image_file, output_image_splits_dir))

