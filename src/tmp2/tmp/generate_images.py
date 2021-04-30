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
import glob

## Variables
## TODO: store in file or statistic varibles, which can share with multiple scripts
SAMPLE      = ['SAMPLE','sample_ID','ID']
CNV_CHR     = ['chr', 'CHR', 'CHROMOSOME', 'chromosome','Chr']
CNV_START   = ['cnv_start', 'start', 'PRED_START', 'START','Start']
CNV_END     = ['cnv_stop', 'stop', 'PRED_END', 'END','End']
CNV_INTERVAL= ['INTERVAL','CNV_INTERVAL']
CNV_TYPE    = ['cnv_type','type','TYPE','CNV', 'CNV_TYPE']
NUM_TARGETS = ['NUM_TARGETS','targets']
GSD_LABEL   = ['LABEL_VAL','label','LABEL','GSD_info']
PRE_LABEL   = ['CNLearn_PRED_LABEL', 'PRED_LABEL']

fixed_win_num = 3 # the number of targets per group
#color_del = (0,1,0) #green
#color_dup = (1,0,0) #red
color_del, color_dup = (0,0,1), (0,0,1) #blue for three classes labels
## Input
RD_norm_dir     = sys.argv[1]
ref_samples_dir = sys.argv[2]
cnv_file        = sys.argv[3] 
output_path     = sys.argv[4]
corr_threshold  = float(sys.argv[5])
try:
    sge_task_id = int(sys.argv[6])
except:
    sge_task_id = 'all'

## Output files and folders
output_EntireCNV_image_dir  = output_path + '/entire_cnvs/'
output_SplitCNV_image_dir   = output_path + '/split_cnvs/'
os.makedirs(output_EntireCNV_image_dir,  exist_ok=True)
os.makedirs(output_SplitCNV_image_dir,   exist_ok=True)

## CNV info
cnv_data_df = pd.read_table(cnv_file, header=0, sep=r'\,|\t', engine='python')

## Functions
def fetchRDdata_byTabix(RD_norm_dir, sampleID, cnv_chr, cnv_start, cnv_end, fixed_win_num):
    # tabix RD file to fetch 
    #RD_norm_file = RD_norm_dir+sampleID+'.cov.bed.norm.gz'
    RD_norm_file = fetch_relative_file_path(RD_norm_dir, sampleID,'gz')
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
    RD_fetched_df_tmp.loc[:, 'fixed_win_num'] = [val for val in np.arange(1,math.ceil((len(RD_fetched_df_tmp)/fixed_win_num))+1) for i in range(fixed_win_num)][0:len(RD_fetched_df_tmp)]
    RD_fetched_df = RD_fetched_df_tmp
    del RD_fetched_df_tmp
    # change the type of columns
    RD_fetched_df[["start"]] = RD_fetched_df[["start"]].astype(int)
    RD_fetched_df[["end"]] = RD_fetched_df[["end"]].astype(int)
    RD_fetched_df[["GC"]] = RD_fetched_df[["GC"]].astype(float)
    RD_fetched_df[["RD_raw"]] = RD_fetched_df[["RD_raw"]].astype(float)
    RD_fetched_df[["RD_norm"]] = RD_fetched_df[["RD_norm"]].astype(float)
    return RD_fetched_df

def loadRefSamplesID(ref_samples_file, corr_threshold):
    ref_samplesID_df = pd.read_table(ref_samples_file,
                                        #low_memory=False,
                                        header=None, sep=' |\t', engine='python',
                                        names=['sampleID', 'r2'])
    # filter by r^2
    result_df = ref_samplesID_df[ref_samplesID_df['r2']>=corr_threshold]
    return result_df

def fetchRefRDdata_byTabix(ref_samples_file, cnv_chr, cnv_start, cnv_end, fixed_win_num, corr_threshold):
    # load reference sample ID
    ref_samplesID_df = loadRefSamplesID(ref_samples_file, corr_threshold)
    # load RD normalized data and fetch RD given the cnv region for each reference sample
    reference_RD_df = pd.DataFrame(columns=['chr', 'start', 'end', 'GC', 'RD_raw', 'RD_norm', 'sample'])
    for index, row in ref_samplesID_df.iterrows():  
        ref_sampleID = row[0]
        try:
            ref_RD_cnv_region = fetchRDdata_byTabix(RD_norm_dir, ref_sampleID, cnv_chr, cnv_start, cnv_end, fixed_win_num)
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

def fetch_relative_file_path(RD_norm_dir, sampleID, suffix):
    sample_rd_file = None
    sample_rd_likely_file = RD_norm_dir+sampleID+'.*.'+suffix
    sample_rd_file_list = glob.glob(sample_rd_likely_file)

    if len(sample_rd_file_list) == 1:
        sample_rd_file = sample_rd_file_list[0]
    else:
        print("   -[Error]: there are multiple files for %s "%(sample_rd_likely_file))
        #TODO: support manually selection of one file
        pdb.set_trace()

    if not os.path.exists(sample_rd_file):
        print("    -[Error]: error in normalized RD file of %s "%(sample_rd_file))
        #TODO: write to log file
        exit(0)
    return sample_rd_file

# ## Processing
## Parse header
cnv_data_header = cnv_data_df.columns.tolist()
col_sampleID  = cnv_data_header.index(fetch_colName(cnv_data_header,SAMPLE))
try:
    col_cnv_interval = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_INTERVAL))
except:
    col_cnv_chr   = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_CHR))
    col_cnv_start = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_START))
    col_cnv_end   = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_END))

col_cnv_type  = cnv_data_header.index(fetch_colName(cnv_data_header,CNV_TYPE))
try:
    col_cnv_num_targets = cnv_data_header.index(fetch_colName(cnv_data_header,['NUM_TARGETS','NUM_TARG']))
except:
    col_cnv_num_targets = None
    pass

try:
    col_cnv_canoes= cnv_data_header.index(fetch_colName(cnv_data_header,['CANOES','CANOES_RT']))
except:
    col_cnv_canoes = None
    pass

try:
    col_cnv_xhmm  = cnv_data_header.index(fetch_colName(cnv_data_header,['XHMM','XHMM_RT']))
except:
    col_cnv_xhmm = None
    pass

try:
    col_cnv_clamms= cnv_data_header.index(fetch_colName(cnv_data_header,['CLAMMS','CLAMMS_RT']))
except:
    col_cnv_clamms = None
    pass

try:
    col_cnv_numCarriers = cnv_data_header.index(fetch_colName(cnv_data_header,['Num_Carriers(inGivenCohort)']))
except:
    col_cnv_numCarriers = None
    pass

try:
    col_GSD_label = cnv_data_header.index(fetch_colName(cnv_data_header,GSD_LABEL))
except:
    col_GSD_label = None
    pass

try:
    col_PRE_label = cnv_data_header.index(fetch_colName(cnv_data_header, PRE_LABEL))
except:
    col_PRE_label = None
    pass

#for index, row in cnv_data_df.iterrows(): 
#     row = next(cnv_data_df.iterrows())[1]
#     index = 1
#    if index + 1 < start_point:
#        continue
index = sge_task_id-1
row   = cnv_data_df.iloc[index]

sampleID  = row[col_sampleID]
try:
    cnv_interval = row[col_cnv_interval]
    cnv_chr, cnv_start, cnv_end = df.parseInterval(cnv_interval) 
except:
    cnv_chr   = row[col_cnv_chr]
    cnv_start = int(row[col_cnv_start])
    cnv_end   = int(row[col_cnv_end])
cnv_type  = row[col_cnv_type]

if cnv_type == 1:
    cnv_type = "DEL"
elif cnv_type == 3:
    cnv_type = "DUP"
else:
    pass

case_sample_color = color_del if cnv_type == 'DEL' else color_dup

try:
    cnv_num_targets = row[col_cnv_num_targets]
except:
    cnv_num_targets =  'NA'
    
try:
    cnv_canoes = str(row[col_cnv_canoes])
except:
    cnv_canoes = 'NA'

try:    
    cnv_xhmm = str(row[col_cnv_xhmm])
except:
    cnv_xhmm = 'NA'

try:
    cnv_clamms = str(row[col_cnv_clamms])
except:
    cnv_clamms = 'NA'

try:
    cnv_num_carriers = str(row[col_cnv_numCarriers])
except:
    cnv_num_carriers = 'NA'

try:    
    cnv_gsd_label = row[col_GSD_label]
except:
    cnv_gsd_label = 'NA'

try:    
    cnv_CNLearn_label = row[col_PRE_label]
except:
    cnv_CNLearn_label = 'NA'

if cnv_gsd_label == 'NA':
    cnv_gsd_str = 'NA'
else:    
    if cnv_gsd_label == 1:
        cnv_gsd_str = "True"
    elif cnv_gsd_label == 0:
        cnv_gsd_str = "False"
    else:
        cnv_gsd_str = cnv_gsd_label 

if cnv_CNLearn_label == 'NA':
    cnv_CNLearn_str = 'NA'
else:    
    cnv_CNLearn_str = "True" if cnv_CNLearn_label == 1 else "False"

print("[%d|%d] Illustrating: %s %s:%d-%d %s #targets:%s Label:%s"% \
       (len(cnv_data_df), index+1, sampleID, cnv_chr, cnv_start, cnv_end, cnv_type, str(cnv_num_targets), cnv_gsd_str))

## Import RD data info
print("  --Step1. Fetching RD data for case sample ...")
RD_cnv_region_df = fetchRDdata_byTabix(RD_norm_dir, sampleID, cnv_chr, cnv_start, cnv_end, fixed_win_num)

## Fetch Read depth data for reference samples in terms of CNV boundary       
print("  --Step2. Fetching RD data for reference samples ...")
#ref_samples_file = ref_samples_dir+sampleID+'.ref.samples.txt.bz2'
ref_samples_file = fetch_relative_file_path(ref_samples_dir, sampleID,'txt')

if not os.path.exists(ref_samples_file):
    print("    -[Error]: error in reference samples related file for %s in %s"%(sampleID, ref_samples_dir))
    exit(0)
reference_RD_df = fetchRefRDdata_byTabix(ref_samples_file, cnv_chr, cnv_start, cnv_end, fixed_win_num, corr_threshold)
    
## plot whole cnv
print("  --Step3. Illustrating an image for the whole CNV ...")
title_info = sampleID+" "+str(cnv_chr)+":"+str(cnv_start)+"-"+str(cnv_end)+" "+cnv_type + \
             " "+ str((cnv_end-cnv_start)/1000) + 'kb'+ " #targets:"+str(cnv_num_targets) + \
             " #wins:" + str(len(RD_cnv_region_df)) + "\nCANOES:"+cnv_canoes + " XHMM:"+ \
             cnv_xhmm + " CLAMMS:"+cnv_clamms + " #Carriers:"+cnv_num_carriers + \
             " CN-Learn:"+cnv_CNLearn_str + "\nGSD_Label:"+cnv_gsd_str

image_file = str(index+1)+"_"+sampleID+"_"+str(cnv_chr)+"_"+str(cnv_start)+"_"+str(cnv_end) + \
             "_"+str(cnv_num_targets)+"tgs_"+str(len(RD_cnv_region_df))+"wins_"+cnv_type+".png"

fig = plt.figure(dpi=150,figsize=(10, 7)) 
ax_rd = fig.subplots(nrows=1, ncols=1)

### plot reference samples
for sample_reader in reference_RD_df["sample"].unique():
            ref_sample_df = reference_RD_df[reference_RD_df["sample"]==sample_reader]
            ax_rd.plot((ref_sample_df["start"]+ref_sample_df["end"])/2, ref_sample_df["RD_norm"],
                        color='grey', marker='.', linewidth=0.2)

### plot case sample
ax_rd.plot((RD_cnv_region_df["start"]+RD_cnv_region_df["end"])/2, RD_cnv_region_df["RD_norm"], \
            color=case_sample_color , marker='o', linewidth=2)
ax_rd.set_title(title_info)
plt.savefig(output_EntireCNV_image_dir+image_file)
plt.close() 

## plot split CNV for each three targets
print("  --Step4. Illustrating images for the CNV splited by each %d windows ..."%fixed_win_num)
split_cnv_path_list = []

cnv_target_num = len(RD_cnv_region_df)
sub_img_num = 0
for i in range(0, cnv_target_num-1, fixed_win_num-1): #for 3 target fixed window, the pace is 2
    sub_img_num += 1
    if i+2 <= cnv_target_num-1:
        start_index = i
        stop_index  = i+2
    elif i+1 <= cnv_target_num-1:
        start_index = i-1
        stop_index  = i+1
    elif i >= cnv_target_num-1:
        break
    print("\tFor %d targets in total, processing targets: %d-%d"%(cnv_target_num, start_index, stop_index))
    # init image info
    split_win = str(int(len(RD_cnv_region_df)/(fixed_win_num-1)))
    title_split_info = title_info + " Images:" + split_win + "-" + str(sub_img_num) 
    image_split_file = str(index+1)+"_"+sampleID+"_"+str(cnv_chr)+"_"+ \
            str(cnv_start)+"_"+str(cnv_end)+"_"+cnv_gsd_str+"_"+cnv_type + \
            "_"+str(cnv_num_targets)+"tgs_" + str(len(RD_cnv_region_df)) + \
            "wins_splits" + split_win + "_" + str(sub_img_num)+".png"

    fig = plt.figure(dpi=150,figsize=(7, 7)) 
    ax_rd = fig.subplots(nrows=1, ncols=1)
    
    # select data
    RD_cnv_region_split = RD_cnv_region_df.iloc[start_index:stop_index+1]

    if len(np.unique(RD_cnv_region_split["chr"]) )>1:
        print("[Error] many chr?")
        pdb.set_trace()

    split_chr = np.unique(RD_cnv_region_split["chr"])[0]
    split_start_from = np.min(RD_cnv_region_split["start"])
    split_start_to   = np.max(RD_cnv_region_split["start"])
    reference_RD_df_split = reference_RD_df[(reference_RD_df["start"]>=split_start_from) &
                                            (reference_RD_df["start"]<=split_start_to)]
    
    # plot reference samples
    for sample_reader in reference_RD_df_split["sample"].unique():
        ref_sample_df =  reference_RD_df_split[reference_RD_df_split["sample"]==sample_reader]
        ax_rd.plot((ref_sample_df["start"]+ref_sample_df["end"])/2, ref_sample_df["RD_norm"],
                    color = 'grey', marker='.', linewidth=0.2)
    
    # plot case sample
    ax_rd.plot((RD_cnv_region_split["start"]+RD_cnv_region_split["end"])/2, RD_cnv_region_split["RD_norm"],
            color=case_sample_color, marker='o', linewidth=2)
    ax_rd.set_title(title_split_info)
    plt.savefig(output_SplitCNV_image_dir+image_split_file)
    plt.close()
    split_cnv_path_list.append(output_SplitCNV_image_dir+image_split_file)
print("  --[Done]. Images have output to %s and %s."%(output_SplitCNV_image_dir+image_file, output_SplitCNV_image_dir))

## address the image path to cnv_info file
cnv_data_df.loc[index, 'split_cnv_path'] = '\n'.join([each_path for each_path in split_cnv_path_list])
