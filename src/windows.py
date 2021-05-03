import pysam
import function as func
import pdb

def calculate_gc(target_list, ref_file):
    fasta_f = pysam.FastaFile(ref_file)
    GC_percentage = []
    chr = target_list[0][0]
    print('Calculating GC content for targets in chr ' + chr)
    
    for i in range(len(target_list)):
        i_chr   = target_list[i][0]
        i_start = int(target_list[i][1])
        i_end   = int(target_list[i][2])

        if(i_chr != chr):
            chr = i_chr
            print('Calculating GC content for targets in chr ' + i_chr)
            
        r_region = fasta_f.fetch(i_chr, i_start, i_end)
        reg_len  = i_end - i_start + 1
        GC = 0
        for b in range(len(r_region)):
            if str(r_region[b]).upper() == 'G' or str(r_region[b]).upper() == 'C':
                GC += 1
        GC_p = round(GC/reg_len, 3)
        target_list[i].extend([GC_p])

    return target_list

def build_windows(target_file, ref_file, output_file):
    output_list = []
    target_list = func.fileToList_tab(target_file)

    for reader in target_list:
        target_chr    = reader[0]
        target_start  = int(reader[1])
        target_end    = int(reader[2])
        target_length = target_end - target_start
        
        # number of windows
        if target_length < 1000:
            win_num = 1
        else:
            win_num = int(target_length/500)

        # windows    
        if win_num == 1:
            win_start = target_start
            win_end   = target_end
            output_list.append([target_chr, win_start, win_end])
        else:
            win_length = int(target_length/win_num+0.5) 
            for win_i in range(0, win_num):
                win_start = target_start + win_i*win_length
                if win_i+1 != win_num:
                    win_end = target_start + win_i*win_length + win_length
                else:
                    win_end = target_end
                output_list.append([target_chr, win_start, win_end])
    
    output_gc_list = calculate_gc(output_list, ref_file)
    func.output_to_file(output_gc_list, output_file)            
