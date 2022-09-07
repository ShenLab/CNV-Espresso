import function as func
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pdb
import os
import math
import pysam
import sys
import random

## Variables
global_var_dict = func.global_variables()
SAMPLE          = global_var_dict['SAMPLE']
CNV_CHR         = global_var_dict['CNV_CHR']
CNV_START       = global_var_dict['CNV_START']
CNV_END         = global_var_dict['CNV_END']
CNV_INTERVAL    = global_var_dict['CNV_INTERVAL']
CNV_TYPE        = global_var_dict['CNV_TYPE']
NUM_TARGETS     = global_var_dict['NUM_TARGETS']
GSD_LABEL       = global_var_dict['CNV_LABEL']
REF_LABEL       = global_var_dict['REF']

def collect_images(cnv_file, output_path, transmission_analysis):
    '''file and folder'''
    if transmission_analysis == True:
        output_EntireCNV_image_dir  = output_path + '/images_transmission_analysis/'
    else:
        output_EntireCNV_image_dir  = output_path + '/images/' 

    if not os.path.exists(cnv_file):
        print("[Error] The cnv list file is not exists. Please check %s"%cnv_file)
        exit(0)
 
    if not os.path.exists(output_EntireCNV_image_dir):
        print("[Error] The image folder is not exists. Please check %s"%output_EntireCNV_image_dir)
        exit(0)
    else:    
        if os.path.splitext(cnv_file)[-1][1:] == 'csv':
            cnv_data_df = pd.read_csv(cnv_file)
        else:    
            cnv_data_df = pd.read_table(cnv_file)

    '''creat columns or fillna'''
    if 'num_of_win' in cnv_data_df.columns:
        cnv_data_df['num_of_win'] = cnv_data_df['num_of_win'].fillna("-")
    else:
        cnv_data_df.insert(cnv_data_df.shape[1], 'num_of_win', "-")
    if 'img_path' in cnv_data_df.columns:
        cnv_data_df['img_path'] = cnv_data_df['img_path'].fillna("-")
    else:
        cnv_data_df.insert(cnv_data_df.shape[1], 'img_path', "-")
    
    if transmission_analysis == True and 'Offspring_img_path' in cnv_data_df.columns:
        cnv_data_df['Offspring_img_path'] = cnv_data_df['Offspring_img_path'].fillna("-")
    if transmission_analysis == True and 'Offspring_img_path' not in cnv_data_df.columns:
        cnv_data_df.insert(cnv_data_df.shape[1], 'Offspring_img_path', "-")

    ## Parse header
    cnv_data_header  = cnv_data_df.columns.tolist()
    col_sampleID     = func.fetch_colName(SAMPLE, cnv_data_header)[0]
    col_cnv_interval = func.fetch_colName(CNV_INTERVAL, cnv_data_header)[0]
    col_cnv_chr      = func.fetch_colName(CNV_CHR, cnv_data_header)[0]
    col_cnv_start    = func.fetch_colName(CNV_START, cnv_data_header)[0]
    col_cnv_end      = func.fetch_colName(CNV_END, cnv_data_header)[0]
    col_cnv_type     = func.fetch_colName(CNV_TYPE, cnv_data_header)[0]
    col_GSD_label    = func.fetch_colName(GSD_LABEL, cnv_data_header)[0]
    col_PRE_label    = func.fetch_colName(['PRE','PRED_LABEL'], cnv_data_header)[0]
    col_cnv_canoes   = func.fetch_colName(['CANOES','CANOES_RT'], cnv_data_header)[0]
    col_cnv_xhmm     = func.fetch_colName(['XHMM','XHMM_RT'], cnv_data_header)[0]
    col_cnv_clamms   = func.fetch_colName(['CLAMMS','CLAMMS_RT'], cnv_data_header)[0]
    col_cnv_num_targets = func.fetch_colName(NUM_TARGETS, cnv_data_header)[0]
    col_cnv_numCarriers = func.fetch_colName(['Num_Carriers(inGivenCohort)'], cnv_data_header)[0]
    col_OffspringID     = func.fetch_colName(['OffspringID'], cnv_data_header)[0]


    # fetch the CNV info from cnv_data_df
    num_cnv = len(cnv_data_df)
    for index, row in cnv_data_df.iterrows():
        sampleID  = row[col_sampleID]
        try:
            cnv_interval = row[col_cnv_interval]
            cnv_chr, cnv_start, cnv_end = func.parseInterval(cnv_interval) 
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

        # if index < 1768:
        #     continue
        ## check if the image is generated/exist and update the file location.
        sample_img_file = func.fetch_relative_file_path(output_EntireCNV_image_dir, 
                                                       "*"+sampleID+"*"+cnv_chr+'*'+str(cnv_start)+"_"+str(cnv_end)+"*"+cnv_type,'png')

        if transmission_analysis == True:
            offspringID        = row[col_OffspringID] if col_OffspringID != None else 'NA'
            offspring_img_file = func.fetch_relative_file_path(output_EntireCNV_image_dir, 
                                                           "*"+offspringID+"*"+cnv_chr+'*'+str(cnv_start)+"_"+str(cnv_end)+"*"+cnv_type,'png')      

        if os.path.exists(str(sample_img_file)):
            print("[%d|%d] Image file exists %s %s:%d-%d %s add the img_path to dataframe."% \
                   (num_cnv, index+1, sampleID, cnv_chr, cnv_start, cnv_end, cnv_type))
            cnv_data_df.loc[index, 'img_path'] = sample_img_file

            if transmission_analysis == True and os.path.exists(str(offspring_img_file)):
                   print("[%d|%d] Image file exists %s %s:%d-%d %s add the offspring img_path to dataframe."% \
                           (num_cnv, index+1, offspringID, cnv_chr, cnv_start, cnv_end, cnv_type))
                   cnv_data_df.loc[index, 'Offspring_img_path'] = offspring_img_file

    cnv_data_df.to_csv(cnv_file, index=False)
    print("Done. output to %s"%cnv_file)
