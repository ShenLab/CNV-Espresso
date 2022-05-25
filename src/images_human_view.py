import function as func
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
import pdb
import os
import re
import math
import vcf
import pysam
import sys
import lockfile
import datetime
from time import sleep
from os.path import exists

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
PRE_LABEL       = global_var_dict['REF']

fixed_win_num = 3 # the number of targets per group
color_del, color_dup = (0,0,1), (0,0,1) #blue for three classes labels

def getBafInfo(snv_vcf_file, sampleID, cnv_chr, figure_left_coordinate, figure_right_coordinate):
    vcf_reader = None
    vcf_reader = vcf.Reader(filename=snv_vcf_file)
    baf_x, baf_y, baf_het, baf_DP = [], [], [], []
    baf_het, baf_hom = [], []

    if cnv_chr not in vcf_reader.contigs:
       if "chr" + cnv_chr in vcf_reader.contigs:
           cnv_chr = "chr" + cnv_chr
       else:
           print("[Error]: there is no chromosome %s in %s"%(str(cnv_chr), snv_vcf_file))
           pdb.set_trace()
    
    if vcf_reader != None:
        records = vcf_reader.fetch(cnv_chr, figure_left_coordinate, figure_right_coordinate + 1)
        snv_num = 0
        record_num = 0
        for record in records:
            record_num += 1
            #if record.is_snp and record.genotype(sampleID).is_variant:
            if record.genotype(sampleID).is_variant:
                snv_num += 1
                snp_pos = record.POS
                snp_het = "r" if record.genotype(sampleID).is_het == True else "b"
                try:
                    snp_DP = record.genotype(sampleID).data.DP
                    ref_allele_reads = record.genotype(sampleID).data.AD[0]
                    alt_allele_reads = sum(x for x in record.genotype(sampleID).data.AD[1: len(record.genotype(sampleID).data.AD)])
                    baf = round((alt_allele_reads+1)/(ref_allele_reads+alt_allele_reads+1),2)
                    baf_x.append(snp_pos)
                    baf_y.append(baf)
                    baf_DP.append(snp_DP)
                    if record.genotype(sampleID).is_het == True:
                        baf_het.append([snp_pos, baf, snp_DP])
                    else:
                        baf_hom.append([snp_pos, baf, snp_DP])
                    time_stamp = datetime.datetime.now()
                except:
                    pass
        print("record_num in vcf:",record_num)
    else:
        print(num, sampleID, 'not in', snv_vcf_file)
        pdb.set_trace()
    # list to pd.dataFrame
    baf_het_df = pd.DataFrame(baf_het, columns=['snv_pos','baf','snv_dp'])
    baf_hom_df = pd.DataFrame(baf_hom, columns=['snv_pos','baf','snv_dp'])
    return([baf_het_df, baf_hom_df])

def draw_baf_figure(ax, snv_vcf_file, cnv_sample, cnv_chr, figure_left_coordinate, figure_right_coordinate, info):
    [baf_het_df, baf_hom_df] = getBafInfo(snv_vcf_file, cnv_sample, cnv_chr, figure_left_coordinate, figure_right_coordinate)

    if info == 'Father':
        baf_color = 'b'
        baf_alpha = 0.4
    elif info == 'Mother':
        baf_color = 'g'
        baf_alpha = 0.4
    elif info == 'Offspring':
        baf_color = 'r'
        baf_alpha = 0.8
    else:
        pdb.set_trace()

    scatter1 = ax.scatter(baf_het_df['snv_pos'], baf_het_df['baf'], sizes=baf_het_df['snv_dp'], color=baf_color, alpha=baf_alpha, marker=".", label="Het")
    scatter2 = ax.scatter(baf_hom_df['snv_pos'], baf_hom_df['baf'], sizes=baf_hom_df['snv_dp'], color=baf_color, alpha=baf_alpha, marker=".", label="Hom")
    
    # rugs; red for het, blue for hom
    if info == 'Offspring':
        sns.rugplot(data=baf_het_df, x='snv_pos', color='r', alpha=.5)
        sns.rugplot(data=baf_hom_df, x='snv_pos', color='b', alpha=.5)

def generate_one_image(vcf_file, cnv_data_df, sge_task_id, col_dict, cnv_info_w_img_file, 
                        RD_norm_dir, ref_samples_dir, output_path, suffix, corr_threshold, flanking, 
                        split_img, trio, ped_file, overwrite_img):
    index = sge_task_id-1
    row   = cnv_data_df.iloc[index]

    # ignore img if overwrite is turned off
    if (overwrite_img == 'False' and cnv_data_df.loc[index, 'img_path'] != "-" and exists(cnv_data_df.loc[index, 'img_path'])):
        print("Image %d existed. Ignore this CNV, since `overwrite_img` is turned off."%index)
        return None

    col_sampleID       = col_dict['col_sampleID']
    col_cnv_interval   = col_dict['col_cnv_interval']
    col_cnv_chr        = col_dict['col_cnv_chr']
    col_cnv_start      = col_dict['col_cnv_start']
    col_cnv_end        = col_dict['col_cnv_end']
    col_cnv_type       = col_dict['col_cnv_type']
    col_cnv_num_targets= col_dict['col_cnv_num_targets']
    col_cnv_num_wins   = col_dict['col_cnv_num_wins']
    col_cnv_canoes     = col_dict['col_cnv_canoes']
    col_cnv_xhmm       = col_dict['col_cnv_xhmm']
    col_cnv_clamms     = col_dict['col_cnv_clamms']
    col_cnv_numCarriers= col_dict['col_cnv_numCarriers']
    col_GSD_label      = col_dict['col_GSD_label']
    col_PRE_label      = col_dict['col_PRE_label']

    output_EntireCNV_image_dir  = output_path + '/images/' 
    output_SplitCNV_image_dir   = output_path + '/images_split_cnvs/'

    sampleID  = row[col_sampleID]
    try:
        cnv_interval = row[col_cnv_interval]
        cnv_chr, cnv_start, cnv_end = func.parseInterval(cnv_interval) 
    except:
        cnv_chr   = row[col_cnv_chr]
        cnv_start = int(row[col_cnv_start])
        cnv_end   = int(row[col_cnv_end])
    cnv_type  = row[col_cnv_type]
    cnv_type  = 'DEL' if (cnv_type == 1 or cnv_type == 'DEL') else 'DUP'
    cnv_type  = re.sub(r'[^\w\s]','',cnv_type)

    if cnv_type == 1:
        cnv_type = "DEL"
    elif cnv_type == 3:
        cnv_type = "DUP"
    else:
        pass
    #case_sample_color = color_del if cnv_type == 'DEL' else color_dup
    case_sample_color = 'red'

    cnv_num_targets   = row[col_cnv_num_targets] if col_cnv_num_targets != None else 'NA'
    cnv_num_wins      = row[col_cnv_num_wins] if col_cnv_num_wins != None else 'NA'
    cnv_canoes        = str(row[col_cnv_canoes]) if col_cnv_canoes != None else 'NA'
    cnv_xhmm          = str(row[col_cnv_xhmm])   if col_cnv_xhmm   != None else 'NA'
    cnv_clamms        = str(row[col_cnv_clamms]) if col_cnv_clamms != None else 'NA'
    cnv_num_carriers  = str(row[col_cnv_numCarriers]) if col_cnv_numCarriers != None else 'NA'
    cnv_gsd_label     = row[col_GSD_label] if col_GSD_label != None else 'NA'
    cnv_CNLearn_label = row[col_PRE_label] if col_PRE_label != None else 'NA'

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
    if str(flanking).upper() == 'TRUE':
        cnv_length   = cnv_end   - cnv_start +1
        figure_left  = max(0, cnv_start - cnv_length)
        figure_right = cnv_end   + cnv_length
    else:
        figure_left  = cnv_start
        figure_right = cnv_end

    print("Generating CNV %s:%d-%d %s with boundary located at %s:%d-%d"%(cnv_chr, cnv_start, cnv_end, cnv_type, 
                                                                        cnv_chr, figure_left, figure_right))

    if str(trio).upper() == 'TRUE':
        try:
            paternalID, maternalID = func.getParentsID(sampleID, ped_file=ped_file)
        except:
            print("In the trio mode, a pedigree file is required. please check.")
            paternalID = None
            maternalID = None
            pdb.set_trace()
            
    ## Import RD and BAF data info
    print("  --Step1. Fetching RD data for case sample ...")
    RD_cnv_region_df = func.fetchRDdata_byTabix(RD_norm_dir, sampleID, cnv_chr, figure_left, figure_right, fixed_win_num, colname='RD_norm')
    if str(trio).upper() == 'TRUE':
        paternal_rd_df = func.fetchRDdata_byTabix(RD_norm_dir, paternalID, cnv_chr, figure_left, figure_right, fixed_win_num, colname='PaternalNormRD')
        maternal_rd_df = func.fetchRDdata_byTabix(RD_norm_dir, maternalID, cnv_chr, figure_left, figure_right, fixed_win_num, colname='MaternalNormRD')

    ## Fetch RD data for reference samples based on CNV/figure boundaries       
    print("  --Step2. Fetching RD data for reference samples ...")
    ref_samples_file = func.fetch_relative_file_path(ref_samples_dir, sampleID,'txt*')
    
    if not os.path.exists(ref_samples_file):
        print("    -[Error]: error in reference samples related file for %s in %s"%(sampleID, ref_samples_dir))
        exit(0)
    reference_RD_df = func.fetchRefRDdata_byTabix(RD_norm_dir, ref_samples_file, 
                                                    cnv_chr, figure_left, figure_right, 
                                                    fixed_win_num, corr_threshold)

    ################################################################
    ## plot the entire cnv into one image
    ################################################################
    print("  --Step3. Illustrating an image for the entire CNV ...")
    RD_cnv_only_region_df = func.fetchRDdata_byTabix(RD_norm_dir, sampleID, cnv_chr, cnv_start, cnv_end, fixed_win_num, colname='RD_norm')
    # Using `RD_cnv_only_region_df.shape[0]` instead of `cnv_num_wins`, since it works for WGS data
    title_info = sampleID+" "+str(cnv_chr)+":"+str(cnv_start)+"-"+str(cnv_end)+" "+cnv_type + \
                 " "+ str((cnv_end-cnv_start)/1000) + 'kb' + " #Targets(Wins):"+ str(RD_cnv_only_region_df.shape[0]) #str(cnv_num_wins)
    
    if suffix==None:
        suffix_content=""
    else:
        suffix_content="_"+str(suffix)

    if str(trio).upper() == 'TRUE':
        image_file = str(index+1).zfill(len(str(cnv_data_df.shape[0])))+"_"+sampleID + \
                     "_"+str(cnv_chr)+"_"+str(cnv_start)+"_"+str(cnv_end)+"_"+cnv_type+suffix_content+"_trio.pdf"
    else:
        image_file = str(index+1).zfill(len(str(cnv_data_df.shape[0])))+"_"+sampleID + \
                     "_"+str(cnv_chr)+"_"+str(cnv_start)+"_"+str(cnv_end)+"_"+cnv_type+suffix_content+".pdf"

    ### Calculate the means and sigmas of ref samples for each target region
    '''
    Note: Uniques are returned in order of appearance. This does NOT sort.
          Here, we assume that there are no two CNVs within a sample have the exactly same start points.
    '''
    target_start_uni = reference_RD_df['start'].unique() 
    target_stop_uni  = reference_RD_df['end'].unique() 
    if len(target_start_uni) != len(target_stop_uni):
        print("[Error] the number of unique `target_start` is not equal to the number of unique `target_stop`. Must be something wrong!")
        pdb.set_trace()
    target_chr_uni   = reference_RD_df['chr'][0:len(target_start_uni)]
    vector = np.vectorize(np.float)
    target_pos_uni   = vector(target_start_uni+ (target_stop_uni - target_start_uni + 1)/2)
    
    ref_ribbon_list = []
    for target_i in range(0,len(target_start_uni)):
        target_i_mu    = reference_RD_df[reference_RD_df['start']==target_start_uni[target_i]]['RD_norm'].mean()
        target_i_sigma = reference_RD_df[reference_RD_df['start']==target_start_uni[target_i]]['RD_norm'].std()
        target_chr   = target_chr_uni[target_i]
        target_start = target_start_uni[target_i]
        target_stop  = target_stop_uni[target_i]
        target_mid   = target_start + (target_stop-target_start+1)/2
        ref_ribbon_list.append([target_chr, target_start, target_stop, target_mid, target_i_mu, target_i_sigma])
    ref_ribbon_df = pd.DataFrame(ref_ribbon_list, columns =['chr', 'start', 'end', 'position', 'mu', 'sigma'])
    
    ### assembly the essential data (target, mu, sigma)
    if str(trio).upper() == 'TRUE':
        case_w_ref_df = pd.merge(ref_ribbon_df, RD_cnv_region_df, on=["chr", "start","end"])
        case_w_ref_df = pd.merge(case_w_ref_df, paternal_rd_df[['chr','start','end','PaternalNormRD']], on=["chr", "start","end"])
        case_w_ref_df = pd.merge(case_w_ref_df, maternal_rd_df[['chr','start','end','MaternalNormRD']], on=["chr", "start","end"])
    else:
        case_w_ref_df = pd.merge(ref_ribbon_df, RD_cnv_region_df)

    ### remove outliers (where mu or sigma equals 0)
    essential_df = case_w_ref_df[(case_w_ref_df['mu']!=0) & (case_w_ref_df['sigma']!=0)]
    #TODO: take care about the sigma outliers
    essential_df = essential_df[(essential_df['RD_norm']-essential_df['mu'])/essential_df['mu'] < 1] 

    ### Preprocessing
    fig   = plt.figure(dpi=150,figsize=(10, 10)) 
    ax_rd, ax_baf = fig.subplots(nrows=2, ncols=1, sharex=True, sharey=False)
    ax_rd.axvspan(xmin=cnv_start, xmax=cnv_end, facecolor='y', alpha=0.1)
    Y_MAX =  1.1
    Y_MIN = -1.1

    ### plot the ribbion
    ax_rd.fill_between(np.array(essential_df['position']), -2*np.array(essential_df['sigma']), 2*np.array(essential_df['sigma']), color='lightgrey', alpha=0.8)
    ax_rd.fill_between(np.array(essential_df['position']), -1*np.array(essential_df['sigma']), 1*np.array(essential_df['sigma']), color='darkgrey',  alpha=0.8)

    ### plot associations
    ax_rd.axhline(y= 0.5, linestyle='dashed', color='dimgrey', alpha=.5)
    ax_rd.axhline(y=-0.5, linestyle='dashed', color='dimgrey', alpha=.5)
    ax_rd.set_ylim([Y_MIN, Y_MAX])
    ax_rd.set_xlim([figure_left, figure_right])
    sns.rugplot(data=essential_df, x="position", color='black', ax=ax_rd)
    ax_rd.set_title(title_info, fontweight='heavy')

    ### plot case sample
    # divide into CNV and boundary regions
    cnv_region_df = essential_df[(essential_df['start']>=cnv_start) & (essential_df['end']<=cnv_end)]
    boundary_region_df = essential_df[(essential_df['start']<cnv_start) | (essential_df['end']>cnv_end)]

    if str(trio).upper() == 'TRUE':
        # illustrate cnv regions (link dots w/ line)
        ax_rd.plot(cnv_region_df['position'], (cnv_region_df['PaternalNormRD']-cnv_region_df['mu'])/cnv_region_df['mu'], \
                    color='b', markerfacecolor='none', marker='.', markersize=10, linewidth=1, alpha=0.4, label='Father')
        ax_rd.plot(cnv_region_df['position'], (cnv_region_df['MaternalNormRD']-cnv_region_df['mu'])/cnv_region_df['mu'], \
                    color='g', markerfacecolor='none', marker='.', markersize=10, linewidth=1, alpha=0.4, label='Mother')
        ax_rd.plot(cnv_region_df['position'], (cnv_region_df['RD_norm']-cnv_region_df['mu'])/cnv_region_df['mu'], \
                    color='r', markerfacecolor='none', marker='.', markersize=10, linewidth=1, alpha=0.8, label='Offspring')
        # illustrate boundary regions (just dots without line)
        ax_rd.plot(boundary_region_df['position'], (boundary_region_df['PaternalNormRD']-boundary_region_df['mu'])/boundary_region_df['mu'], \
                    color='b', markerfacecolor='none', marker='.', markersize=10, linewidth=1, linestyle = 'None', alpha=0.25, label='Father')
        ax_rd.plot(boundary_region_df['position'], (boundary_region_df['MaternalNormRD']-boundary_region_df['mu'])/boundary_region_df['mu'], \
                    color='g', markerfacecolor='none', marker='.', markersize=10, linewidth=1, linestyle = 'None', alpha=0.25, label='Mother')
        ax_rd.plot(boundary_region_df['position'], (boundary_region_df['RD_norm']-boundary_region_df['mu'])/boundary_region_df['mu'], \
                    color='r', markerfacecolor='none', marker='.', markersize=10, linewidth=1, linestyle = 'None', alpha=0.4, label='Offspring')
        ax_rd.legend(loc='upper right')

    else:
        # illustrate cnv regions (link dots w/ line)
        ax_rd.plot(cnv_region_df['position'], (cnv_region_df['RD_norm']-cnv_region_df['mu'])/cnv_region_df['mu'], \
                    color='r', markerfacecolor='none', marker='.', markersize=10, linewidth=1, alpha=0.8)
    
        # illustrate boundary regions (just dots without line)
        ax_rd.plot(boundary_region_df['position'], (boundary_region_df['RD_norm']-boundary_region_df['mu'])/boundary_region_df['mu'], \
                    color='r', markerfacecolor='none', marker='.', markersize=10, linewidth=1, linestyle = 'None', alpha=0.4)
    
    ax_rd.set_ylabel('Norm depth relative to diploid', fontsize=15, fontweight='normal')
    ax_rd.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, symbol='%', is_latex=False))

    ## BAF plot
    ### plot the line/area in background
    ax_baf.axhline(y=0.5, color='dimgrey', linestyle='solid',  linewidth=1, alpha=0.5)
    ax_baf.axhline(y=2/3, color='dimgrey', linestyle='dashed', linewidth=1, alpha=0.3)
    ax_baf.axhline(y=1/3, color='dimgrey', linestyle='dashed', linewidth=1, alpha=0.3)
    ax_baf.axvspan(xmin=cnv_start, xmax=cnv_end, facecolor='y', alpha=0.1)

    ### plot the SNVs
    boundary_length = round((figure_right - figure_left+1)/1000, 2)
    print("          Getting SNV information and drawing BAF plot. Boundary for BAF: %s:%d-%d, %fkb ..."%(cnv_chr, figure_left, figure_right, boundary_length))
    draw_baf_figure(ax_baf, vcf_file, sampleID, cnv_chr, figure_left, figure_right, info='Offspring')
    if str(trio).upper() == 'TRUE':
        draw_baf_figure(ax_baf, vcf_file, paternalID, cnv_chr, figure_left, figure_right, info='Father')
        draw_baf_figure(ax_baf, vcf_file, maternalID, cnv_chr, figure_left, figure_right, info='Mother')
    ax_baf.set_xlim(figure_left, figure_right)
    ax_baf.set_ylim(-0.05, 1.05)
    ax_baf.set_xlabel('Chromosome position', fontsize=15, fontweight='normal')
    ax_baf.set_ylabel('Alternative allele fraction', fontsize=15, fontweight='normal')
    ax_baf.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, symbol='%', is_latex=False))

    fig.align_labels() 
    ### write the img path to the cnv_file_w_img_file
    print("  --Step4. Output image file to %s."%(output_EntireCNV_image_dir+image_file))
    img_path = output_EntireCNV_image_dir+image_file
    plt.savefig(img_path)
    plt.close()  

    print("  --Step5. Update the %s with img path."%cnv_info_w_img_file)
    lock_flag = lockfile.LockFile(cnv_info_w_img_file)
    while lock_flag.is_locked():
        sleep(0.1)
    lock_flag.acquire()
    cnv_data_df = pd.read_csv(cnv_info_w_img_file)
    cnv_data_df.loc[index, 'img_path'] = img_path
    cnv_data_df.to_csv(cnv_info_w_img_file, index=False)
    lock_flag.release()
    print("  --[Done].")

        
def generate_images_human_view(RD_norm_dir, ref_samples_dir, cnv_file, vcf_file, output_path, suffix, corr_threshold,  
                               flanking, split_img, sge_task_id, trio, pedigree_file, overwrite_img):
    try:
        sge_task_id = int(sge_task_id)
    except:
        sge_task_id = 'all'

    ## Output images and folders
    output_EntireCNV_image_dir  = output_path + '/images/' 
    os.makedirs(output_EntireCNV_image_dir, exist_ok=True)
    if split_img == True:
        output_SplitCNV_image_dir   = output_path + '/images_split_cnvs/'
        os.makedirs(output_SplitCNV_image_dir, exist_ok=True)
    
    ## Output cnv info file with image path
    '''
    Idea: at the very beginning, just copy the cnv_file to cnv_w_img_file, then
          once the cnv_w_img_file existed, insert/update each img path into the corresponding cell.
          Note: to avoid the file writing conflict, we input and output it w/ img path at the end.
    '''
    ## CNV info
    cnv_info_w_img_file = output_path + '/cnv_info_w_img.csv'
    if not os.path.exists(cnv_info_w_img_file):
        file_dialect = func.fetchFileDialect(cnv_file)
        if file_dialect == '\t':
            cnv_data_df =  pd.read_table(cnv_file,low_memory=False,header=0, sep='\t')
        elif file_dialect == ',':
            cnv_data_df = pd.read_csv(cnv_file)
        cnv_data_df.to_csv(cnv_info_w_img_file, index=False)
    else:
        cnv_data_df = pd.read_csv(cnv_info_w_img_file)

    ## Parse header
    cnv_data_header  = cnv_data_df.columns.tolist()
    col_sampleID     = func.fetch_colName(SAMPLE, cnv_data_header)[0]
    col_cnv_interval = func.fetch_colName(CNV_INTERVAL, cnv_data_header)[0]
    col_cnv_chr      = func.fetch_colName(CNV_CHR, cnv_data_header)[0]
    col_cnv_start    = func.fetch_colName(CNV_START, cnv_data_header)[0]
    col_cnv_end      = func.fetch_colName(CNV_END, cnv_data_header)[0]
    col_cnv_type     = func.fetch_colName(CNV_TYPE, cnv_data_header)[0]
    col_GSD_label    = func.fetch_colName(GSD_LABEL, cnv_data_header)[0]
    col_PRE_label    = func.fetch_colName(PRE_LABEL, cnv_data_header)[0]
    col_cnv_canoes   = func.fetch_colName(['CANOES','CANOES_RT'], cnv_data_header)[0]
    col_cnv_xhmm     = func.fetch_colName(['XHMM','XHMM_RT'], cnv_data_header)[0]
    col_cnv_clamms   = func.fetch_colName(['CLAMMS','CLAMMS_RT'], cnv_data_header)[0]
    col_cnv_num_targets = func.fetch_colName(NUM_TARGETS, cnv_data_header)[0]
    col_cnv_num_wins = func.fetch_colName(['num_of_win'], cnv_data_header)[0]
    col_cnv_numCarriers = func.fetch_colName(['Num_Carriers(inGivenCohort)'], cnv_data_header)[0]

    col_dict = {
        'col_sampleID'     : col_sampleID, 
        'col_cnv_interval' : col_cnv_interval,
        'col_cnv_chr'      : col_cnv_chr,
        'col_cnv_start'    : col_cnv_start,
        'col_cnv_end'      : col_cnv_end,
        'col_cnv_type'     : col_cnv_type,
        'col_cnv_num_targets': col_cnv_num_targets,
        'col_cnv_num_wins' : col_cnv_num_wins,
        'col_cnv_canoes'   : col_cnv_canoes,
        'col_cnv_xhmm'     : col_cnv_xhmm,
        'col_cnv_clamms'   : col_cnv_clamms,
        'col_cnv_numCarriers': col_cnv_numCarriers,
        'col_GSD_label'    : col_GSD_label,
        'col_PRE_label'    : col_PRE_label
    }
    if sge_task_id == False:
        for index, row in cnv_data_df.iterrows(): 
           index += 1
           generate_one_image(vcf_file, cnv_data_df, index, col_dict, cnv_info_w_img_file, 
                            RD_norm_dir, ref_samples_dir, output_path, suffix, corr_threshold, flanking, 
                            split_img, trio, pedigree_file, overwrite_img)

    else:
        generate_one_image(vcf_file, cnv_data_df, sge_task_id, col_dict, cnv_info_w_img_file,
                            RD_norm_dir, ref_samples_dir, output_path, suffix, corr_threshold, flanking, 
                            split_img, trio, pedigree_file, overwrite_img)

