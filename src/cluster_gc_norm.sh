#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=100:0:0
#$ -l h_vmem=18G
#S -l mem_free=18G
#$ -pe smp 1

#SBATCH -o job.%j.out
#SBATCH -J RD_normalization
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 18G
#SBATCH --time=99:99:99

set -oe pipefail

echo "Job started on `hostname` at `date`"

ID=${SGE_TASK_ID}${SLURM_ARRAY_TASK_ID}
SCRIPT_DIR=$1
WINDOWS_FILE=$2
INPUT_RD_FILE_LIST_W_PATH=$3
OUTPUT_DIR=$4

if [ ! -d $OUTPUT_DIR ]
then
    mkdir -p $OUTPUT_DIR
fi

RD_file=$(head -n $ID $INPUT_RD_FILE_LIST_W_PATH | tail -n 1 | awk '{print $1}')
file_name=${RD_file%.*}
file_name=$(basename $file_name)
#file_name=$(basename $RD_file .bz2)

echo ${RD_file}
python ${SCRIPT_DIR}/cnv_espresso.py normalization \
    --windows ${WINDOWS_FILE} \
    --input ${RD_file} \
    --output ${OUTPUT_DIR}

#bgzip ${OUTPUT_DIR}${file_name}".norm"
#tabix -p bed ${OUTPUT_DIR}${file_name}".norm.gz"
echo "Job ended on `hostname` at `date`"
