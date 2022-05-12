import os
import sys
import pysam 
import numpy as np
import function as func
import pdb

def gc_normalization(windows_file, input_file_list, output_dir, debug_flag):
    print('Loading windows.bed ...')
    #TODO: columns of loadWindows function
    windows_dict = func.loadWindows(windows_file)
    windows_chr  = windows_dict['chr']
    windows_start= windows_dict['start']
    windows_stop = windows_dict['stop']
    windwos_gc   = windows_dict['gc']
    #windows_mappability = windows_dict['mappability']

    for input_file_reader in input_file_list:
        input_file = input_file_reader[0]
        (input_dir,input_name) = os.path.split(input_file)
        (input_filename,input_extension) = os.path.splitext(input_name)

        output_file = output_dir + input_filename + '.norm'

        print('Loading %s ...'% input_file)
        input_sample_dict = func.loadRD(input_file)
        input_sample_chr   = input_sample_dict['chr']
        input_sample_start = input_sample_dict['start']
        input_sample_stop  = input_sample_dict['stop']
        input_sample_rd    = input_sample_dict['RD']

        if not np.array_equal(input_sample_chr, windows_chr) or \
                not np.array_equal(input_sample_start, windows_start) or \
                not np.array_equal(input_sample_stop, windows_stop):
                    print("[Error] The windows file is not consisitant with input sample RD file.")
                    pdb.set_trace()
                    sys.exit(0)

        GC_percentage = np.round(windwos_gc*100)
        GC_index = {}
        for ind in range(len(GC_percentage)):
            gc = int(GC_percentage[ind])

            if gc in GC_index:
                GC_index[gc].append(ind)
            else:
                GC_index[gc] = [ind]

        print('Normalizing by GC percentage...')
        corrected_rd   = np.zeros(len(input_sample_rd), dtype=float)
        overall_median = np.median(GC_percentage)
        overall_mean   = np.mean(GC_percentage)

        sum_num_targets = 0
        for gc in GC_index.keys():
            t_ind = GC_index[gc]
            t_median = np.median(input_sample_rd[t_ind])
            sum_num_targets += len(t_ind)
            if debug_flag == True:
                print("GC:",gc,"t-median:",t_median,"#target:",len(t_ind))

            if t_median == 0:
                if debug_flag == True:
                    print('[WARNING] Median read depth of [GC=%d] targets is 0. Set %d RD of these target to 0.'%(len(t_ind),gc))
                corrected_rd[t_ind] = 0
            else:
                corrected_rd[t_ind] = (input_sample_rd[t_ind] * overall_median / t_median) / overall_mean

        if debug_flag == True:
            print("For across all %d targets, overall_median: %f, overall_mean:%f"%(sum_num_targets, overall_median,  overall_mean))
        corrected_rd = np.around(corrected_rd,decimals=6)

        # Output
        output_ndarray = np.transpose(np.array([input_sample_chr, input_sample_start, input_sample_stop, \
                GC_percentage, input_sample_rd, corrected_rd]))
        print('Saving and tabix normalized read depth to file %s.gz'%output_file)
        func.output_to_file(output_ndarray, output_file)                
        pysam.tabix_index(output_file, preset="bed",force=True)

