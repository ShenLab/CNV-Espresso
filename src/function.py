from __future__ import division
import vcf, os
import csv
import re
import pdb
import datetime
import numpy as np
import pandas as pd
import pysam
import sys
import glob
import math

def main():
    pass
# input/output/parse files
def global_variables():
    global_variable_dict = {
        "SAMPLE"      : ['SAMPLE','SAMPLE_ID','ID','sample_name','sampleID'],
        "CNV_INTERVAL": ['CNV_INTERVAL','INTERVAL'],
        "CNV_CHR"     : ['CHR', 'CNV_CHR', 'CHROMOSOME', 'CHROMOSOMES'],
        "CNV_START"   : ['CNV_START', 'PRED_START', 'START','st_bp'],
        "CNV_END"     : ['CNV_STOP', 'STOP', 'PRED_END', 'END', 'ed_bp'],
        "CNV_TYPE"    : ['CNV_TYPE','TYPE','CNV'],
        "NUM_TARGETS" : ['NUM_TARG','NUM_TARGETS','Num_Targets_Wins','num_target'],
        "CNV_LABEL"   : ['LABEL_VAL','LABEL','GSD_info'],
        "REF"         : ['ref','REF','Ref','Reference','hg19','hg38']
    }
    return global_variable_dict

def norm_cnv_type(cnv_type):
    if cnv_type.upper() in ["DEL", "DELETION", "RARE_DEL", "RARE_DELETION", "RARE DELETION"]:
        cnv_type = "DEL"
    elif cnv_type.upper() in ["DUP", "DUPLICATION", "RARE_DUP", "RARE_DUPLICATION", "RARE DUPLICATION"]:
        cnv_type = "DUP"
    else:
        pass

    return cnv_type

def parseInterval(cnv_interval):
    pattern = re.compile(r'(\S+):(\d+)-(\d+)')
    try:
        loc_list = pattern.match(cnv_interval)
        cnv_chr = loc_list.group(1)
        cnv_start = int(loc_list.group(2))
        cnv_stop = int(loc_list.group(3))
    except:
        print("[ERROR] interval does not match 'chr:start-stop' pattern.")
        cnv_chr = None
        cnv_start = None
        cnv_stop = None
    return [cnv_chr, cnv_start, cnv_stop]

def output_to_file(results,file_name, show_time=False):
    path = os.path.dirname(file_name)
    if not os.path.exists(path) and path!='':
        os.mkdir(path)

    fp = open(file_name,'w')
    if type(results) == list:
        for result_line in results:
            result_line = [str(x) for x in result_line]
            fp.write("\t".join(result_line) + os.linesep)            

    elif type(results) == dict:
        for key in results.keys():
            for result_line in results[str(key)]:
                result_line = [str(x) for x in result_line]
                fp.write("\t".join(result_line) + os.linesep)
    elif type(results) == np.ndarray:
        np.savetxt(file_name, results, fmt='%s', delimiter="\t")
    else:
        print('[ERROR]: Unsupport results type!')
        return
    
    print('[INFO]: File outputs to %s'%file_name)
    if show_time != False:
        time_stamp = datetime.datetime.now()
        print(time_stamp.strftime('[%Y.%m.%d-%H:%M:%S]'))
    fp.close()

def fileToList(file_name):
    result_list = [] 
    extension = os.path.splitext(file_name)[-1][1:]
    if extension == 'bz2':
        fp = bz2.open(file_name, "rt")
    else:
        fp = open(file_name)
    for row in fp:
        row = row.strip()
        row = row.replace("\"","")
        row = re.split(r'[;,\t]', row)
        result_list.append(row)
    fp.close()  
    return result_list 

def fileToList_tab(file_name):
    result_list = [] 
    fp = open(file_name)
    for row in fp:
        row = row.strip()
        row = row.replace("\"","")
        row = re.split('\t',row)
        result_list.append(row)
    fp.close()  
    return result_list 

def fetch_colName(keyWord_list, colName_list):
    keyWord_list = [col_name.upper() for col_name in keyWord_list]
    #colName_list = [col_name.upper() for col_name in colName_list]

    index_result = []
    keyWord_result = []
    for index, keyWord in enumerate(colName_list):
        if keyWord.upper() in keyWord_list:
            index_result.append(index)
            #keyWord_result.append(keyWord)
            keyWord_result.append(colName_list[index])
    if len(keyWord_result) == 1:
        return index_result[0], keyWord_result[0]
    elif len(keyWord_result) > 1:
        print("[Error] Multiple keywords for one. Please take a look")
        pdb.set_trace()
    else:    
        return [None, None]

def parse_pedigreeFile(file_name):
    unaffected_parents_list = []
    try:
        fp = open(file_name,'r')
        row_num = 0
        ped_indiv_dict = {}
        ped_family_dict = {}
        for reader in fp:
            row_num += 1
            reader = reader.strip()
            reader = reader.split('\t')
            FamilyID = reader[0]
            IndividualID = reader[1]
            PaternalID = reader[2]
            MaternalID = reader[3]
            Sex = reader[4]
            Phenotype = reader[5]
            if len(reader) == 7:
                Role = reader[6]
                ped_indiv_dict[IndividualID] = [FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role]
                if FamilyID in ped_family_dict:
                    ped_family_dict[FamilyID].append([FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role])
                else:
                    ped_family_dict[FamilyID] = [[FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role]]
            elif len(reader) >= 8:
                Role = reader[6]
                Ancestry = reader[7]
                ped_indiv_dict[IndividualID] = [FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role, Ancestry]
                if FamilyID in ped_family_dict:
                    ped_family_dict[FamilyID].append([FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role, Ancestry])
                else:
                    ped_family_dict[FamilyID] = [[FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype, Role, Ancestry]]
            else:
                ped_indiv_dict[IndividualID] = [FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype]
                if FamilyID in ped_family_dict:
                    ped_family_dict[FamilyID].append([FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype])
                else:
                    ped_family_dict[FamilyID] = [[FamilyID, IndividualID, PaternalID, MaternalID, Sex, Phenotype]]
        time_stamp = datetime.datetime.now()
        print(time_stamp.strftime('%Y.%m.%d-%H:%M:%S'),"\tPedigree file %s imported."%file_name)
        return [ped_indiv_dict, ped_family_dict]
    except:
        print('Error df.parse_pedigreeFile function when openning %s'%file_name)
        return [None, None]

def getParentsID(sampleID, ped_file):
    [ped_indiv_dict, ped_family_dict] = parse_pedigreeFile(ped_file)
    paternalID = ped_indiv_dict[sampleID][2]
    maternalID = ped_indiv_dict[sampleID][3]
    return([paternalID, maternalID])

def getHeader(input_file):
    header_list = []
    fp = open(input_file)
    for reader in fp:
        reader = reader.strip()
        reader = reader.replace('"','')
        reader = re.split('\t| ',reader)
        header_list = reader
        break
    fp.close()
    return header_list

def loadRD(filename):
    rd_matrix = np.loadtxt(filename, delimiter='\t', skiprows=0,
            dtype={'names': ('chr', 'start', 'stop','RD'),
                'formats': ('U5', np.int, np.int, np.float)}) 
    rd_matrix_sort = np.sort(rd_matrix, order=['chr', 'start'] )
    return rd_matrix_sort

def loadNormRD(filename):
    if os.path.splitext(filename)[-1][1:] == "gz":
        f = gzip.open(filename,mode='rt')
    elif os.path.splitext(filename)[-1][1:] == "bz2":
        f = bz2.open(filename,mode='rt')
    else:
        f = open(filename)
    chr_str = np.loadtxt(f, dtype=np.str, delimiter='\t', skiprows=0, usecols=(0,)) 
    f.seek(0)
    start = np.loadtxt(f, dtype=np.int, delimiter='\t', skiprows=0, usecols=(1,)) 
    f.seek(0)
    stop = np.loadtxt(f, dtype=np.int, delimiter='\t', skiprows=0, usecols=(2,)) 
    f.seek(0)
    GC = np.loadtxt(f, dtype=np.float, delimiter='\t', skiprows=0, usecols=(3,)) 
    f.seek(0)
    RD_raw = np.loadtxt(f, dtype=np.float, delimiter='\t', skiprows=0, usecols=(4,))
    f.seek(0)
    RD_norm = np.loadtxt(f, dtype=np.float, delimiter='\t', skiprows=0, usecols=(5,))
    return {'chr': chr_str, 'start': start, 'stop': stop, 'GC': GC, 'RD_raw': RD_raw, 'RD_norm': RD_norm}

def loadWindows(windows_filename):
    try:
        targets = np.loadtxt(windows_filename, delimiter='\t', skiprows=0,
                dtype={'names': ('chr', 'start', 'stop','interval','unknow1','gc','mappability','unknow2'),
                       'formats': ('U5', np.int, np.int, 'U20', np.int, np.float, np.float, np.int)})
    except:
        try:
            targets = np.loadtxt(windows_filename, delimiter='\t', skiprows=0,
                     dtype={'names': ('chr', 'start', 'stop','gc'),
                         'formats': ('U5', np.int, np.int,  np.float)})
        except:
            print("[Error] Please check the windows.bed file.")
            pdb.set_trace()
    targets_sort = np.sort(targets, order=['chr', 'start'])        
    return targets_sort

def fetch_norm_rd(sampleID, sample_norm_file):
    #df = pd.read_table(sample_norm_file,low_memory=False,header=None, sep='\t',
    #                   names=['chr', 'start','end','GC','RD',sampleID])
    df = pd.read_table(sample_norm_file,low_memory=False,header=None, sep='\t')
    if df.shape[1] == 6: # CNV-Espresso's format
        df.columns = ['chr', 'start','end','GC','RD',sampleID]
    elif df.shape[1] == 4: # CLAMMS' format
        df.columns = ['chr', 'start','end',sampleID]
    else:
        pdb.set_trace()  
    return df[sampleID].to_frame()

def fetch_sampleID_from_filepath(filepath):
    filepath,tempfilename = os.path.split(filepath[0])
    #sampleID = tempfilename.replace(".cov.bed.norm.gz","")
    sampleID = tempfilename.replace("coverage","")
    sampleID = sampleID.replace("cov","")
    sampleID = sampleID.replace("bed","")
    sampleID = sampleID.replace("norm","")
    sampleID = sampleID.replace("gz","")
    sampleID = sampleID.replace(".","")
    return sampleID

def fetchRDdata_byTabix(RD_norm_dir, sampleID, cnv_chr, cnv_start, cnv_end, fixed_win_num, colname):
    '''
    Aim: Fetch RD signal from tabix file 
    [Update] 2021.05.18 compatible with both CNV-Espresso and CLAMMS's normlized file
    [Update] 2021.07.07 add colname as a variable
    '''
    RD_norm_file = fetch_relative_file_path(RD_norm_dir, sampleID,'gz')
    if not os.path.exists(RD_norm_dir):
        print('No tabular file: %s'%RD_norm_file) 
        return ["No tabular file"]
    if RD_norm_file == None:
        print("No this RD_norm_file %s for sample:%s"%(str(RD_norm_file), sampleID))    
        return ["No norm file"]

    if not os.path.exists(RD_norm_file+'.tbi'):
        pysam.tabix_index(RD_norm_file, seq_col=0, start_col=1, end_col=2) # Need to add '-p bed'
    # fetch
    f = pysam.TabixFile(RD_norm_file)
    try:
        RD_fetched_data = f.fetch(cnv_chr, int(cnv_start), int(cnv_end), parser=pysam.asTuple())
    except:
        try:
            cnv_chr = "chr"+cnv_chr
            RD_fetched_data = f.fetch(cnv_chr, int(cnv_start), int(cnv_end), parser=pysam.asTuple())
        except:
            print("[Error] Failed to fetch SNVs in %s:%d-%d for sample:%s"%(str(cnv_chr), int(cnv_start), int(cnv_end), sampleID))
            pdb.set_trace()

    RD_fetched_df = pd.DataFrame(data=RD_fetched_data)
    if RD_fetched_df.shape[1] == 6: # for CNV-Espresso normlized file
        RD_fetched_df.columns = ['chr', 'start', 'end', 'GC', 'RD_raw', colname]
    elif RD_fetched_df.shape[1] == 4: # for CLAMMS normlized file
        RD_fetched_df.columns = ['chr', 'start', 'end', colname]
    elif RD_fetched_df.shape[1] == 0:
        print("    --[Warnning] No target in %s %s:%d-%d"%(sampleID, cnv_chr, int(cnv_start), int(cnv_end)))
        return RD_fetched_df
    else:
        print("    --[Error] in the format of input file.") 
        pdb.set_trace()

    # add a new column as target groups
    RD_fetched_df_tmp = RD_fetched_df.copy()
    RD_fetched_df_tmp.loc[:, 'fixed_win_num'] = [val for val in np.arange(1,math.ceil((len(RD_fetched_df_tmp)/fixed_win_num))+1) for i in range(fixed_win_num)][0:len(RD_fetched_df_tmp)]
    RD_fetched_df = RD_fetched_df_tmp
    del RD_fetched_df_tmp

    # change the type of columns
    RD_fetched_df[["start"]]   = RD_fetched_df[["start"]].astype(int)
    RD_fetched_df[["end"]]     = RD_fetched_df[["end"]].astype(int)
    RD_fetched_df[[colname]]   = RD_fetched_df[[colname]].astype(float)
    
    try:
        RD_fetched_df[["GC"]]     = RD_fetched_df[["GC"]].astype(float)
        RD_fetched_df[["RD_raw"]] = RD_fetched_df[["RD_raw"]].astype(float)
    except:
        # Norm files from CLAMMS do not contain 'GC' and 'RD_raw' columns
        pass
    return RD_fetched_df

def fetchFileDialect(input_file):
    header_list = []
    extension = os.path.splitext(input_file)[-1][1:]
    if extension == 'bz2':
        fp = bz2.open(input_file, "rt")
    else:
        fp = open(input_file)

    dialect = csv.Sniffer().sniff(fp.readline())
    fp.close()
    return(dialect.delimiter)

def loadRefSamplesID(ref_samples_file, corr_threshold):
    ref_samplesID_df = pd.read_table(ref_samples_file,
                                        #low_memory=False,
                                        header=None, sep=' |\t', engine='python',
                                        names=['sampleID', 'r2'])
    # filter by r^2
    result_df = ref_samplesID_df[ref_samplesID_df['r2']>=corr_threshold]
    return result_df

def fetchRefRDdata_byTabix(RD_norm_dir, ref_samples_file, cnv_chr, cnv_start, cnv_end, fixed_win_num, corr_threshold):
    '''
    [Update] 2021.05.18 compatible with both CNV-Espresso and CLAMMS's normlized file
    '''
    # load reference sample ID
    ref_samplesID_df = loadRefSamplesID(ref_samples_file, corr_threshold)
    # load RD normalized data and fetch RD given the cnv region for each reference sample
    ref_espresso_RD_df = pd.DataFrame(columns=['chr', 'start', 'end', 'GC', 'RD_raw', 'RD_norm', 'sample'])
    ref_clamms_RD_df   = pd.DataFrame(columns=['chr', 'start', 'end', 'RD_norm', 'sample'])

    for index, row in ref_samplesID_df.iterrows():  
        ref_sampleID = row[0]
        try:
            ref_RD_cnv_region = fetchRDdata_byTabix(RD_norm_dir, ref_sampleID, cnv_chr, cnv_start, cnv_end, \
                                                    fixed_win_num,'RD_norm')
            # add a new column as sampleID
            RD_cnv_region_tmp = ref_RD_cnv_region.copy()
            RD_cnv_region_tmp.loc[:, 'sample'] = [ref_sampleID]*len(ref_RD_cnv_region)
            ref_RD_cnv_region = RD_cnv_region_tmp
            del RD_cnv_region_tmp
            # combine results
            if ref_RD_cnv_region.shape[1] == 6+2: # for CNV-Espresso
                ref_espresso_RD_df = ref_espresso_RD_df.append(ref_RD_cnv_region)
            elif ref_RD_cnv_region.shape[1] == 4+2: # for CLAMMS
                ref_clamms_RD_df = ref_clamms_RD_df.append(ref_RD_cnv_region)
            elif ref_RD_cnv_region.shape[1] == 1: # no target in this region
                pass
            else:
                print("[Error] please check the norm file. ")
                pdb.set_trace()
        except:
            print("    -[Error]: error in normalized reference RD file of %s in %s"%(ref_sampleID, RD_norm_dir))
    if ref_espresso_RD_df.shape[0] != 0 and ref_clamms_RD_df.shape[0] == 0:
        # change the type of columns
        ref_espresso_RD_df[["start"]]   = ref_espresso_RD_df[["start"]].astype(int)
        ref_espresso_RD_df[["end"]]     = ref_espresso_RD_df[["end"]].astype(int)
        ref_espresso_RD_df[["RD_norm"]] = ref_espresso_RD_df[["RD_norm"]].astype(float)
        
        try:
            ref_espresso_RD_df[["GC"]]     = ref_espresso_RD_df[["GC"]].astype(float)
            ref_espresso_RD_df[["RD_raw"]] = ref_espresso_RD_df[["RD_raw"]].astype(float)
        except:
            # Norm files from CLAMMS do not contain 'GC' and 'RD_raw' columns
            pass
        return ref_espresso_RD_df
    elif ref_espresso_RD_df.shape[0] == 0 and ref_clamms_RD_df.shape[0] != 0:
        # change the type of columns
        ref_clamms_RD_df[["start"]]   = ref_clamms_RD_df[["start"]].astype(int)
        ref_clamms_RD_df[["end"]]     = ref_clamms_RD_df[["end"]].astype(int)
        ref_clamms_RD_df[["RD_norm"]] = ref_clamms_RD_df[["RD_norm"]].astype(float)
        
        try:
            ref_clamms_RD_df[["GC"]]     = ref_clamms_RD_df[["GC"]].astype(float)
            ref_clamms_RD_df[["RD_raw"]] = ref_clamms_RD_df[["RD_raw"]].astype(float)
        except:
            # Norm files from CLAMMS do not contain 'GC' and 'RD_raw' columns
            pass
        return ref_clamms_RD_df
    elif ref_espresso_RD_df.shape[0] == 0 and ref_clamms_RD_df.shape[0] == 0:
        print("[Error] no reference sample or this is an outlier sample which has extremely low correlation with ref samples?")
        return ref_espresso_RD_df
    else:
        print("[Error] mixed norm files form both CNV-Espresso and CLAMMS ?")
        pdb.set_trace()
    
def fetch_relative_file_path(RD_norm_dir, sampleID, suffix):
    sample_rd_file = None
    #sample_rd_likely_file = RD_norm_dir+'/'+sampleID+'.*.'+suffix
    sample_rd_likely_file = RD_norm_dir+'/'+sampleID+'*.'+suffix
    sample_rd_file_list = glob.glob(sample_rd_likely_file)
    if len(sample_rd_file_list) == 1:
        sample_rd_file = sample_rd_file_list[0]
    else:
        print("   -[Error]: there are zero or multiple files for %s "%(sample_rd_likely_file))

    try:
        if not os.path.exists(sample_rd_file):
            print("    -[Error]: error in normalized RD file of %s "%(sample_rd_file))
            # exit(0)
    except:
        pass
            
    return sample_rd_file

def try_except(fn, default):
    try:
        return fn
    except:
        return default

def showDateTime(end='\n'):
    time_stamp = datetime.datetime.now()
    print(time_stamp.strftime('[%Y.%m.%d-%H:%M:%S]'),end=end)

def extractFilePathNameExtension(full_path):
    path, filename_w_ext = os.path.split(full_path)
    filename, file_extension = os.path.splitext(filename_w_ext)
    return path, filename, file_extension    
