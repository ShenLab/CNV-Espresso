########################################################################
# Script : extract_rc_ratio_cluster.sh                                 #
# Author : Renjie Tan                                                  #
# Date   : 8/25/2020                                                   #
#                                                                      #
# This script calculates read count ratio for each sample by cluster.  #
# The read count was calculated by CANOES algorithm                    # 
#                                                                      # 
########################################################################
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=100:0:0
#$ -l h_vmem=45G
#$ -l mem_free=45G
#$ -o '/share/terra/rt2776/CNV_Espresso/oe/'
#$ -e '/share/terra/rt2776/CNV_Espresso/oe/'

echo "Job started on `hostname` at `date`"
SCRIPTS_DIR=$1
DATA_DIR=$2
SAMPLE_W_BATCH_LIST=$3
OUTPUT_DIR=$4

## Test argu
#echo $0
#echo ${SCRIPTS_DIR}
#echo ${DATA_DIR}
#echo ${SAMPLE_W_BATCH_LIST}
#echo ${OUTPUT_DIR}

TEST_SAMPLE=$(head -n $SGE_TASK_ID $SAMPLE_W_BATCH_LIST | tail -n 1 | awk '{print $1}')
BATCH_NAME=$(head -n $SGE_TASK_ID $SAMPLE_W_BATCH_LIST | tail -n 1 | awk '{print $2}')

if [ ! -f ${OUTPUT_DIR}${TEST_SAMPLE}.txt.gz ];
then
    Rscript --vanilla ${SCRIPTS_DIR}extract_rc_ratio.R ${DATA_DIR} ${TEST_SAMPLE} ${BATCH_NAME} ${OUTPUT_DIR}
    bgzip ${OUTPUT_DIR}${TEST_SAMPLE}.txt
    tabix -p bed ${OUTPUT_DIR}${TEST_SAMPLE}.txt.gz
fi
echo "Job ended on `hostname` at `date`"
