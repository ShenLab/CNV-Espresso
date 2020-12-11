########################################################################
# Script : extract_rc_ratio_sample_check.sh                            #
# Author : Renjie Tan                                                  #
# Date   : 8/26/2020                                                   #
#                                                                      #
# This script checks the extracted read count file for each sample.    #
# Missing samples should be fixed and extracted again properly.        #
########################################################################
#!/bin/bash

echo "Checking process started on `hostname` at `date`"
echo "Missing samples will be listed in below:"
SAMPLE_W_BATCH_LIST=$1
OUTPUT_DIR=$2

num=0
while read LINE  
do  
    let "num++"
    TEST_SAMPLE=$(echo ${LINE} | cut -d " " -f 1)
    if [ ! -f ${OUTPUT_DIR}${TEST_SAMPLE}.txt.gz ];
    then
        echo "[$dt] Missing ${TEST_SAMPLE} on job ${num}. "
    fi
done < ${SAMPLE_W_BATCH_LIST}  

echo "Job ended on `hostname` at `date`"
