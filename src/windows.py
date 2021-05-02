import function as func
import pdb

def build_windows(target_file, ref_file, output_file):
    print(target_file)
    print(ref_file)
    print(output_file)
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
            print(target_chr, target_start, target_end, target_length, "win:", win_start, win_end)
        else:
            win_length = int(target_length/win_num) 
            for win_i in range(0, win_num):
                win_start = target_start + win_i*win_length
                if win_i+1 != win_num:
                    win_end = target_start + win_i*win_length + win_length
                else:
                    win_end = target_end

                print(target_chr, target_start, target_end, target_length, "win:", win_start, win_end)
            pdb.set_trace()
