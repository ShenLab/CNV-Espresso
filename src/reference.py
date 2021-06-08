import function as func
import pandas as pd
import numpy as np
import pdb

'''
TBD: 
        1. Support to take an exist corr_matrix. There is no need to calculate the matrix for each time.
        2. Output zipped ref.sample.txt directly. and the following functions should support to take zipped file directly.
'''

def reference_selection(project_path, norm_list_file, num_ref, corr_threshold):
        func.showDateTime()
        print("[Step1] Parsing the normlized sample list and importing the normlized RD data for each sample ...")
        num = 0
        combined_list = []
        sampleID_list = []
        sample_norm_rd_list = func.fileToList_tab(norm_list_file)

        for sample_rd_path in sample_norm_rd_list:
            sampleID = func.fetch_sampleID_from_filepath(sample_rd_path)
            sample_rd_file = sample_rd_path[0]
            num += 1
            if num % 100 == 0:
                func.showDateTime('\t')
                print("Importing No.%d sample:%s from %s"%(num, sampleID, sample_rd_file))
            rd_df = func.fetch_norm_rd(sampleID, sample_rd_file)
            combined_list.append(rd_df.to_numpy())
            sampleID_list.append(sampleID)

        ## convert to nparray to accelerate the speed
        combined_np_array = np.hstack(combined_list)
        combined_df       = pd.DataFrame(combined_np_array,columns=sampleID_list)
        func.showDateTime()

        # correlation matrix
        func.showDateTime()
        print("[Step2] Calculating the correlation matrix ...")
        corrMatrix = pd.DataFrame(np.corrcoef(combined_np_array, rowvar=False), columns=sampleID_list)

        ## output correlation matrix
        sampleID_str = '\t'.join([sampleID for sampleID in sampleID_list])
        corr_matrix_file = project_path+'/correlation_matrix.gz'
        np.savetxt(corr_matrix_file, corrMatrix, delimiter="\t", header=sampleID_str, comments='')

        ## heatmap
        # sn.heatmap(corrMatrix, annot=True)
        # plt.show()

        func.showDateTime()
        print("[Step3] Selecting reference samples for each case sample ...")
        for case_sampleID in sampleID_list:
            ref_sample_list = []
            ref_sample_df = corrMatrix[case_sampleID].sort_values(ascending=False)
            ref_sample_size = min(num_ref,len(corrMatrix))
            for i in range(1,ref_sample_size+1):
                ref_sampleID = sampleID_list[ref_sample_df.index[i]]
                ref_sample_corr = ref_sample_df.iloc[i]
                if ref_sample_corr >= corr_threshold:
                    ref_sample_list.append([ref_sampleID, ref_sample_corr])

            output_ref_file = project_path + '/ref_samples/'+case_sampleID+'.ref.samples.txt'
            func.output_to_file(ref_sample_list, output_ref_file)
