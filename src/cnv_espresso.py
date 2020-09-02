import os
import re
import pdb
import datetime
import sys
import argparse
import function as df

KEY_WORDS_SAMPLE = ['SAMPLE','sample_ID']
KEY_WORDS_CHR    = ['chr', 'CHR', 'CHROMOSOME', 'chromosome']
KEY_WORDS_START  = ['cnv_start', 'start', 'PRED_START']
KEY_WORDS_END    = ['cnv_stop', 'stop', 'PRED_END']
KEY_WORDS_TYPE   = ['cnv_type','type','TYPE','CNV']

def split_cnv(args):
    input_file  = args.input[0]
    output_file = args.output[0]
    tensor_size = int(args.num_target[0])
    
    print("Import %s ..."%input_file)
    cnv_list = df.fileToList(input_file)
    result_list = []
    result_list.append(cnv_list[0])

    # fetch the column name
    try:
        cnv_sample_col = cnv_list[0].index(df.fetch_colName(KEY_WORDS_SAMPLE, cnv_list[0]))
        if 'INTERVAL' in cnv_list[0]:
            cnv_interval_col = cnv_list[0].index('INTERVAL')
        else:
            cnv_chr_col   = cnv_list[0].index(df.fetch_colName(KEY_WORDS_CHR, cnv_list[0]))
            cnv_start_col = cnv_list[0].index(df.fetch_colName(KEY_WORDS_START, cnv_list[0]))
            cnv_stop_col  = cnv_list[0].index(df.fetch_colName(KEY_WORDS_END, cnv_list[0]))
        
        cnv_type_col    = cnv_list[0].index(df.fetch_colName(KEY_WORDS_TYPE, cnv_list[0]))
        num_targets_col = cnv_list[0].index(df.fetch_colName(['NUM_TARGETS'], cnv_list[0]))
        rc_ratio_col    = cnv_list[0].index(df.fetch_colName(['RC_Ratio'], cnv_list[0]))
    except IOError as e:
        print('[ERROR] Please check your column name.')

    # Main loops
    num=0
    total_num = len(cnv_list[1:])
    for reader in cnv_list[1:]:
        num+=1
        cnv_chr     = reader[cnv_chr_col]
        cnv_start   = int(reader[cnv_start_col])
        cnv_stop    = int(reader[cnv_stop_col])
        cnv_type    = reader[cnv_type_col]
        sampleID    = reader[cnv_sample_col]
        num_targets = int(reader[num_targets_col])
        rc_ratio    = reader[rc_ratio_col]
        
        rc_ratio      = re.sub(r'[\[\]]','',rc_ratio)
        rc_ratio_list = re.split(r'[;,\s]\s*', rc_ratio)
        
        print("---------------------------------------------------")
        num_win = num_targets//tensor_size
        #if sampleID == 'SP0000362' and cnv_start == 103571602:
        #    pdb.set_trace()
        for win_item in range(0, num_win):
            #some of cnv may not have enough rc ratio values. such as SP0000362 1:103571602-103754814
            if tensor_size*win_item+tensor_size <= len(rc_ratio_list): 
                win_rc_list = rc_ratio_list[tensor_size*win_item : tensor_size*win_item+tensor_size]
                result_tmp  = reader[:-1]
                result_tmp.extend([','.join(win_rc_list)])
                result_list.append(result_tmp)
                print(num, num_win, cnv_chr, cnv_start, cnv_stop, cnv_type, sampleID, num_targets, win_rc_list)
            
    df.output_to_file(result_list, file_name = output_file)

def main():
    pass

if __name__ == '__main__':
    main()
    
    VERSION = "0.0.1"
    parser = argparse.ArgumentParser(prog="cnv_espresso", description="This is CNV_Espresso %s, designed to filter CNVs identified from WES data by deep learning." % VERSION)
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    subparsers = parser.add_subparsers(help='Command to be run.')

    # Commands
    cnv_parser = subparsers.add_parser('split', help='Split read count ratios in CNVs into same-size window (default:3 targets)')
    cnv_parser.add_argument('--input',action='store', required=True, metavar='input file',nargs=1,  help="CNV file in bed-like format")
    cnv_parser.add_argument('--output',action='store', required=True, metavar='output file',nargs=1,  help="Output file")
    cnv_parser.add_argument('--num_target',action='store', required=False, default=3, metavar=None, nargs=1, help='Window size')
    cnv_parser.set_defaults(func=split_cnv)
    
    args = parser.parse_args()
    args.func(args)
