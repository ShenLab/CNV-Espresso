import os
import numpy as np
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
import dask.dataframe as dd

import function as func
import itertools
from datetime import datetime 

# # Functions
def fetch_norm_rd(sampleID, sample_norm_file):
    df = pd.read_table(sample_norm_file,low_memory=False,header=None, sep='\t',
                       names=['chr', 'start','end','GC','RD',sampleID])
    return df[sampleID].to_frame()

def fetch_sampleID_from_filepath(filepath):
    filepath,tempfilename = os.path.split(filepath[0])
    sampleID = tempfilename.replace(".cov.bed.norm.gz","")
    return sampleID


# # Variables
num_ref = 100
corr_threshold = -1


# In[5]:


project_path = '/home/rt2776/1000GP/cnv_espresso/'


# In[6]:


sample_norm_rd_file='/home/rt2776/1000GP/cnv_espresso/sample_norm_rd.list'
sample_norm_rd_file='/home/rt2776/1000GP/cnv_espresso/sample_norm_rd_NA12878.list'


# # Processes

# In[7]:


sample_norm_rd_list = func.fileToList_tab(sample_norm_rd_file)


# In[9]:


func.showDateTime()
num = 0
combined_list = []
sampleID_list = []

for sample_rd_path in sample_norm_rd_list:
    sampleID = fetch_sampleID_from_filepath(sample_rd_path)
    sample_rd_file = sample_rd_path[0]
    num += 1
#     print(num, sampleID, sample_rd_file)
    if num % 10 == 0:
        func.showDateTime('\t')
        print("Importing No.%d sample:%s from %s"%(num, sampleID, sample_rd_file))
    rd_df = fetch_norm_rd(sampleID, sample_rd_file)
    combined_list.append(rd_df.to_numpy())
    sampleID_list.append(sampleID)

combined_np_array = np.hstack(combined_list) # convert to nparray to accelerate the speed.
combined_df = pd.DataFrame(combined_np_array,columns=sampleID_list)
# combined_df.columns = sampleID_list

func.showDateTime()


# In[10]:


# faster way than `corrMatrix = combined_df.corr()`
# ref: https://stackoverflow.com/questions/56628363/computing-correlation-matrix-faster-in-pandas
## corrMatrix = pd.DataFrame(np.corrcoef(combined_df.values, rowvar=False), columns=combined_df.columns)

corrMatrix = pd.DataFrame(np.corrcoef(combined_np_array, rowvar=False), columns=sampleID_list)
func.showDateTime()


# In[11]:


corrMatrix


# In[12]:


sampleID_str = '\t'.join([sampleID for sampleID in sampleID_list])


# In[13]:


corr_matrix_file = project_path+'/correlation_matrix.txt'
np.savetxt(corr_matrix_file,corrMatrix, delimiter="\t", header=sampleID_str, comments='')


# In[57]:


# sn.heatmap(corrMatrix, annot=True)
# plt.show()


# # Select reference samples

# In[85]:


for case_sampleID in sampleID_list:
    num_ref=100
    ref_sample_list = []

    print("case sampleID:",case_sampleID)
    ref_sample_df = corrMatrix[case_sampleID].sort_values(ascending=False)
    ref_sample_size = min(num_ref,len(corrMatrix))
    for i in range(1,ref_sample_size):
        ref_sampleID = sampleID_list[ref_sample_df.index[i]]
        ref_sample_corr = ref_sample_df.iloc[i]
        if ref_sample_corr >= corr_threshold:
            ref_sample_list.append([ref_sampleID, ref_sample_corr])

    output_ref_file = project_path + '/ref_samples/'+case_sampleID+'.ref.samples.txt'
    func.output_to_file(ref_sample_list, output_ref_file)


# In[14]:


num_ref=100
ref_sample_list = []
case_sampleID = 'NA12878'

print("case sampleID:",case_sampleID)
ref_sample_df = corrMatrix[case_sampleID].sort_values(ascending=False)
ref_sample_size = min(num_ref,len(corrMatrix))
for i in range(1,ref_sample_size):
    ref_sampleID = sampleID_list[ref_sample_df.index[i]]
    ref_sample_corr = ref_sample_df.iloc[i]
    if ref_sample_corr >= corr_threshold:
        ref_sample_list.append([ref_sampleID, ref_sample_corr])

output_ref_file = project_path + '/ref_samples/'+case_sampleID+'.ref.samples.txt'
func.output_to_file(ref_sample_list, output_ref_file)

