#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=100:0:0
#$ -l h_vmem=18G
#S -l mem_free=18G
#$ -pe smp 1
#$ -o '/share/terra/rt2776/tmp/oe/'
#$ -e '/share/terra/rt2776/tmp/oe/'
echo "Job started on `hostname` at `date`"

ID=$SGE_TASK_ID
WINDOWS_FILE=$1
INPUT_RD_FILE_LIST_W_PATH=$2
OUTPUT_DIR=$3

if [ ! -d $OUTPUT_DIR ]
then
    mkdir -p $OUTPUT_DIR
fi

RD_file=$(head -n $ID $INPUT_RD_FILE_LIST_W_PATH | tail -n 1 | awk '{print $1}')

echo ${RD_file}
python /home/rt2776/cnv_espresso/src/cnv_espresso.py normalization \
    --windows ${WINDOWS_FILE} \
    --input ${RD_file} \
    --output ${OUTPUT_DIR}

echo "Job ended on `hostname` at `date`"
