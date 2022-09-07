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
import lockfile
from time import sleep

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

fixed_win_num = 3 # the number of targets per group
color_del, color_dup = (0,0,1), (0,0,1) #blue for three classes labels

def generate_one_image(cnv_data_df, sge_task_id, col_dict, cnv_info_w_img_file, 
                        RD_norm_dir, ref_samples_dir, output_path, corr_threshold, 
                        flanking, single_img_info, overwrite_img, offspring_img):
    index = int(sge_task_id)-1

    # ignore img if overwrite is turned off
    if (overwrite_img == 'False' and cnv_data_df.loc[index, 'img_path'] != "-"):
        if offspring_img == True:
            if 'Offspring_img_path' in cnv_data_df.columns:
                if cnv_data_df.loc[index, 'Offspring_img_path'] != "-" :
                    print("Ignore %d CNV. Image and corresponding offspring image existed. Meanwhile, `overwrite_img` is turned off."%index)
                    return None
        else:        
            print("Image %d existed. Ignore this CNV, since `overwrite_img` is turned off."%index)
            return None
       
    # header
    row   = cnv_data_df.iloc[index]
    col_sampleID       = col_dict['col_sampleID']
    col_cnv_interval   = col_dict['col_cnv_interval']
    col_cnv_chr        = col_dict['col_cnv_chr']
    col_cnv_start      = col_dict['col_cnv_start']
    col_cnv_end        = col_dict['col_cnv_end']
    col_cnv_type       = col_dict['col_cnv_type']
    col_cnv_num_targets= col_dict['col_cnv_num_targets']
    col_cnv_canoes     = col_dict['col_cnv_canoes']
    col_cnv_xhmm       = col_dict['col_cnv_xhmm']
    col_cnv_clamms     = col_dict['col_cnv_clamms']
    col_cnv_numCarriers= col_dict['col_cnv_numCarriers']
    col_GSD_label      = col_dict['col_GSD_label']
    col_PRE_label      = col_dict['col_PRE_label']
    col_OffspringID    = col_dict['col_OffspringID']

    # output
    output_SplitCNV_image_dir   = output_path + '/images_split_cnvs/'
    if offspring_img == True:
        output_EntireCNV_image_dir  = output_path + '/images_transmission_analysis/'
    else:
        output_EntireCNV_image_dir  = output_path + '/images/'

    if single_img_info == True:
        output_single_img_info_dir  = output_path + '/single_img_info/'

    # fetch the CNV info from cnv_data_df
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

    case_sample_color = color_del if cnv_type == 'DEL' else color_dup

    cnv_num_targets   = row[col_cnv_num_targets] if col_cnv_num_targets != None else 'NA'
    cnv_canoes        = str(row[col_cnv_canoes]) if col_cnv_canoes != None else 'NA'
    cnv_xhmm          = str(row[col_cnv_xhmm])   if col_cnv_xhmm   != None else 'NA'
    cnv_clamms        = str(row[col_cnv_clamms]) if col_cnv_clamms != None else 'NA'
    cnv_num_carriers  = str(row[col_cnv_numCarriers]) if col_cnv_numCarriers != None else 'NA'
    cnv_gsd_label     = row[col_GSD_label] if col_GSD_label != None else 'NA'
    cnv_CNLearn_label = row[col_PRE_label] if col_PRE_label != None else 'NA'
    if offspring_img == True:
        offspringID = row[col_OffspringID] if col_OffspringID != None else 'NA'

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

    ## Confirm the boundries
    if flanking == True or flanking == 'True':
        cnv_length   = cnv_end - cnv_start + 1
        figure_left  = cnv_start - cnv_length/2
        figure_right = cnv_end + cnv_length/2
    else:
        figure_left  = cnv_start
        figure_right = cnv_end

    print(" Generate CNV %s:%d-%d %s with boundary located at %s:%d-%d"%(cnv_chr, cnv_start, cnv_end, cnv_type, 
                                                                        cnv_chr, figure_left, figure_right))
    ## Import RD data info
    print("  --Step1. Fetching RD data for case sample ...")
    RD_cnv_region_df = func.fetchRDdata_byTabix(RD_norm_dir, sampleID, cnv_chr, 
                                                figure_left, figure_right, 
                                                fixed_win_num, colname='RD_norm')
    if offspring_img == True:
        print("  --Step1_offspring. Fetching RD data for the offspring of case sample ...")
        offspring_RD_cnv_region_df = func.fetchRDdata_byTabix(RD_norm_dir, offspringID, cnv_chr, 
                                                figure_left, figure_right, 
                                                fixed_win_num, colname='RD_norm')


    ## Fetch Read depth data for reference samples in terms of CNV boundary       
    print("  --Step2. Fetching RD data for reference samples ...")
    ref_samples_file = func.fetch_relative_file_path(ref_samples_dir, sampleID, '*')
    
    if not os.path.exists(str(ref_samples_file)):
        print("    -[Error]: error in reference samples related file for %s in %s"%(sampleID, ref_samples_dir))
        return None
    reference_RD_df = func.fetchRefRDdata_byTabix(RD_norm_dir, ref_samples_file, 
                                                    cnv_chr, figure_left, figure_right, 
                                                    fixed_win_num, corr_threshold)
    if offspring_img == True:
        print("  --Step2_offspring. Fetching RD data for the offspring reference samples ...")
        offspring_ref_samples_file = func.fetch_relative_file_path(ref_samples_dir, offspringID, '*')
        offspring_reference_RD_df  = func.fetchRefRDdata_byTabix(RD_norm_dir, offspring_ref_samples_file, 
                                                                    cnv_chr, figure_left, figure_right, 
                                                                    fixed_win_num, corr_threshold)

    ## plot the entire cnv into one image
    print("  --Step3. Illustrating an image for the entire CNV ...")
    title_info = sampleID + "  " + str(cnv_chr) + ":" + str(cnv_start)+ "-" + str(cnv_end) + \
                 "  " + str((cnv_end-cnv_start)/1000) + 'kb' + "  #targets:" + str(cnv_num_targets) + \
                 "  #wins:" + str(len(RD_cnv_region_df)) + "  "+ cnv_type

    fig = plt.figure(dpi=150,figsize=(10, 7)) 
    ax_rd = fig.subplots(nrows=1, ncols=1)

    ### plot reference samples
    if not reference_RD_df.empty: 
        for sample_reader in reference_RD_df["sample"].unique():
            ref_sample_df = reference_RD_df[reference_RD_df["sample"]==sample_reader]

            #### Transform data (log differnece(X-axis) and log Y-axis)
            ref_sample_pos_df     = (ref_sample_df["start"]+ref_sample_df["end"])/2
            ref_pos_df_diff       = np.ediff1d(ref_sample_pos_df, to_begin=0)
            ref_pos_df_cumLogDiff = ref_sample_pos_df[0] + np.cumsum(np.log1p(ref_pos_df_diff))
            ax_rd.plot(ref_pos_df_cumLogDiff, np.log1p(ref_sample_df["RD_norm"]), color='grey', marker='.', linewidth=0.2)

    ### plot case sample
    #### Transform data (log differnece(X-axis) and log Y-axis)
    if not RD_cnv_region_df.empty:
        case_sample_pos_df     = (RD_cnv_region_df["start"]+RD_cnv_region_df["end"])/2
        case_pos_df_diff       = np.ediff1d(case_sample_pos_df, to_begin=0)
        case_pos_df_cumLogDiff = case_sample_pos_df[0] + np.cumsum(np.log1p(case_pos_df_diff))
        ax_rd.plot(case_pos_df_cumLogDiff, np.log1p(RD_cnv_region_df["RD_norm"]), color=case_sample_color , marker='o', linewidth=2)

    ax_rd.set_title(title_info)

    ### write the img path to the cnv_file_w_img_file
    image_file = output_EntireCNV_image_dir + str(index+1).zfill(len(str(len(cnv_data_df)))) + "_" + sampleID + "_" + \
                    str(cnv_chr)+"_"+str(cnv_start)+"_"+str(cnv_end) + \
                    "_"+str(cnv_num_targets)+"tgs_"+str(len(RD_cnv_region_df))+"wins_"+cnv_type+".png"                
    image_file = image_file.replace("//", "/")

    if single_img_info == True:
        single_cnv_info_file = output_single_img_info_dir + func.extractFilePathNameExtension(image_file)[1]+'.csv'
        single_cnv_info_file = single_cnv_info_file.replace("//","/")

    print("  --Step4. Output image file to %s."%(image_file))
    plt.savefig(image_file)
    plt.close()  

    '''
    illustrate img for offspring of the case sample
    '''
    if offspring_img == True:
        print("  --Step3_offspring. Illustrating an image for the entire CNV of offspring ...")
        offspring_title_info = offspringID + "  " + str(cnv_chr) + ":" + str(cnv_start)+ "-" + str(cnv_end) + \
                                   "  " + str((cnv_end-cnv_start)/1000) + 'kb' + "  #targets:" + str(cnv_num_targets) + \
                                    "  #wins:" + str(len(RD_cnv_region_df)) + "  "+ cnv_type

        offspring_fig   = plt.figure(dpi=150,figsize=(10, 7)) 
        offspring_ax_rd = offspring_fig.subplots(nrows=1, ncols=1)

        ### plot reference samples
        if not offspring_reference_RD_df.empty: 
            for sample_reader in offspring_reference_RD_df["sample"].unique():
                offspring_ref_sample_df = offspring_reference_RD_df[offspring_reference_RD_df["sample"]==sample_reader]

                #### Transform data (log differnece(X-axis) and log Y-axis)
                offspring_ref_sample_pos_df     = (offspring_ref_sample_df["start"]+offspring_ref_sample_df["end"])/2
                offspring_ref_pos_df_diff       = np.ediff1d(offspring_ref_sample_pos_df, to_begin=0)
                offspring_ref_pos_df_cumLogDiff = offspring_ref_sample_pos_df[0] + np.cumsum(np.log1p(offspring_ref_pos_df_diff))
                offspring_ax_rd.plot(offspring_ref_pos_df_cumLogDiff, np.log1p(offspring_ref_sample_df["RD_norm"]), color='grey', marker='.', linewidth=0.2)

        ### plot case sample
        #### Transform data (log differnece(X-axis) and log Y-axis)
        if not offspring_RD_cnv_region_df.empty:
            offspring_case_sample_pos_df     = (offspring_RD_cnv_region_df["start"]+offspring_RD_cnv_region_df["end"])/2
            offspring_case_pos_df_diff       = np.ediff1d(offspring_case_sample_pos_df, to_begin=0)
            offspring_case_pos_df_cumLogDiff = offspring_case_sample_pos_df[0] + np.cumsum(np.log1p(offspring_case_pos_df_diff))
            offspring_ax_rd.plot(offspring_case_pos_df_cumLogDiff, np.log1p(offspring_RD_cnv_region_df["RD_norm"]), color=case_sample_color , marker='o', linewidth=2)

        offspring_ax_rd.set_title(offspring_title_info)
        ### write the img path to the cnv_file_w_img_file_transmission
        offspring_image_file = output_EntireCNV_image_dir + str(index+1).zfill(len(str(len(cnv_data_df)))) + "_" + offspringID + "_" + \
                                str(cnv_chr)+"_"+str(cnv_start)+"_"+str(cnv_end) + \
                                "_"+str(cnv_num_targets)+"tgs_"+str(len(RD_cnv_region_df))+"wins_"+cnv_type+".png"
        offspring_img_path = offspring_img_path.replace("//","/")
        print("  --Step4_offspring. Output corrsponding offspring image file to %s."%(offspring_img_path))
        plt.savefig(offspring_img_path)
        plt.close()  

    print("  --Step5. Output CNV info with img path: %s"%(cnv_info_w_img_file if single_img_info==False else single_cnv_info_file))
    cnv_data_df.loc[index, 'num_of_win'] = len(RD_cnv_region_df)
    cnv_data_df.loc[index, 'img_path']   = image_file
    if offspring_img == True:
       cnv_data_df.loc[index,'Offspring_img_path'] = offspring_img_path
    
    if single_img_info == False:
        cnv_data_df.to_csv(cnv_info_w_img_file, index=False)
    else:
        if index == 0:
            cnv_data_df.loc[index].to_frame().T.to_csv(single_cnv_info_file, index=False, header=True)
        else:
            cnv_data_df.loc[index].to_frame().T.to_csv(single_cnv_info_file, index=False, header=False)

    print("  --[Done].")

def generate_images(RD_norm_dir, ref_samples_dir, cnv_file, output_path, corr_threshold, flanking, single_img_info, 
                    sge_task_id, job_start, overwrite_img, offspring_img):
    print("sge_task_id:", sge_task_id)
    try:
        sge_task_id = int(sge_task_id)
    except:
        sge_task_id = 'all'

    ## Prepare for output images and folders
    if offspring_img == True:
        output_EntireCNV_image_dir  = output_path + '/images_transmission_analysis/'
    else:
        output_EntireCNV_image_dir  = output_path + '/images/' 
    os.makedirs(output_EntireCNV_image_dir, exist_ok=True)

    ## Prepare for cnv info file with image path
    '''
    Idea: for the very beginning, copy the cnv_file to cnv_w_img_file, then
          once the cnv_w_img_file exists, insert/update each img path into the corresponding cell.
          Note: to avoid the file writing conflict, we input and output it w/ img path at the end.
    '''
    if offspring_img == True:
        cnv_info_w_img_file = output_path + '/cnv_info_w_img_for_tranmission_analysis.csv'
    else:
        cnv_info_w_img_file = output_path + '/cnv_info_w_img.csv'
        
    if not os.path.exists(cnv_info_w_img_file):
        if os.path.splitext(cnv_file)[-1][1:] == 'csv':
            cnv_data_df = pd.read_csv(cnv_file)
        else:    
            cnv_data_df = pd.read_table(cnv_file)
        cnv_data_df.to_csv(cnv_info_w_img_file, index=False)
    else:
        cnv_data_df = pd.read_csv(cnv_info_w_img_file)  

    '''when we run under cluster, we will generate split info files w/ single cnv info with img path.'''
    if single_img_info == True:
        output_single_img_info_dir = output_path + '/single_img_info/'
        os.makedirs(output_single_img_info_dir, exist_ok=True)

    '''creat columns or fillna'''
    if 'num_of_win' in cnv_data_df.columns:
        cnv_data_df['num_of_win'] = cnv_data_df['num_of_win'].fillna("-")
    else:
        cnv_data_df.insert(cnv_data_df.shape[1], 'num_of_win', "-")
    if 'img_path' in cnv_data_df.columns:
        cnv_data_df['img_path'] = cnv_data_df['img_path'].fillna("-")
    else:
        cnv_data_df.insert(cnv_data_df.shape[1], 'img_path', "-")
    
    if offspring_img == True and 'Offspring_img_path' in cnv_data_df.columns:
        cnv_data_df['Offspring_img_path'] = cnv_data_df['Offspring_img_path'].fillna("-")
    if offspring_img == True and 'Offspring_img_path' not in cnv_data_df.columns:
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

    col_dict = {
        'col_sampleID'     : col_sampleID,
        'col_cnv_interval' : col_cnv_interval,
        'col_cnv_chr'      : col_cnv_chr,
        'col_cnv_start'    : col_cnv_start,
        'col_cnv_end'      : col_cnv_end,
        'col_cnv_type'     : col_cnv_type,
        'col_cnv_num_targets': col_cnv_num_targets,
        'col_cnv_canoes'   : col_cnv_canoes,
        'col_cnv_xhmm'     : col_cnv_xhmm,
        'col_cnv_clamms'   : col_cnv_clamms,
        'col_cnv_numCarriers': col_cnv_numCarriers,
        'col_GSD_label'    : col_GSD_label,
        'col_PRE_label'    : col_PRE_label,
        'col_OffspringID'  : col_OffspringID
    }

    if sge_task_id == False:
        if job_start == False:
            for index, row in cnv_data_df.iterrows(): 
               index += 1
               generate_one_image(cnv_data_df, index, col_dict, cnv_info_w_img_file, 
                                    RD_norm_dir, ref_samples_dir, output_path, corr_threshold, 
                                    flanking, single_img_info, overwrite_img, offspring_img)
        else:
            for index in range(int(job_start), len(cnv_data_df)+1):
                generate_one_image(cnv_data_df, index, col_dict, cnv_info_w_img_file, 
                                    RD_norm_dir, ref_samples_dir, output_path, corr_threshold,
                                    flanking, single_img_info, overwrite_img, offspring_img)

                print("  There are %d images left. Continue ..."%(len(cnv_data_df)-index))
    else:
        generate_one_image(cnv_data_df, sge_task_id, col_dict, cnv_info_w_img_file,
                            RD_norm_dir, ref_samples_dir, output_path, corr_threshold,
                             flanking, single_img_info, overwrite_img, offspring_img)

