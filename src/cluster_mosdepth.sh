#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=100:0:0
#$ -l h_vmem=4G
#S -l mem_free=4G
#$ -pe smp 1
set -oe pipefail

#SBATCH -o job.%j.out
#SBATCH -J mosdepth
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBTACH --mem 4G
#SBATCH --time=99:99:99

echo "Job started on `hostname` at `date`"

ID=${SGE_TASK_ID}${SLURM_ARRAY_TASK_ID}
CRAM_FILE_LIST_W_PATH=$1
REF_GENOME=$2
TARGET_PROBES=$3
OUTPUT_RD_DIR=$4

echo "Your input cram files with path: ${CRAM_FILE_LIST_W_PATH}"
echo "Number of processing samples: ${ID}"
echo "Your input Ref_genome file: ${REF_GENOME}"
echo "Your input target probes file: ${TARGET_PROBES}"
echo "Output RD folder: ${OUTPUT_RD_DIR}"

if [ ! -d $OUTPUT_RD_DIR ]
then
    mkdir -p $OUTPUT_RD_DIR
fi

cram_file=$(head -n $ID $CRAM_FILE_LIST_W_PATH | tail -n 1 |awk '{print $1}')
sample_name=$(head -n $ID $CRAM_FILE_LIST_W_PATH | tail -n 1 |awk '{print $2}')
mosdepth -n --fasta ${REF_GENOME} --by ${TARGET_PROBES} --mapq 30 ${OUTPUT_RD_DIR}/${sample_name} $cram_file 

ncol=$(head -n1 <(zcat ${OUTPUT_RD_DIR}/${sample_name}.regions.bed.gz) | awk '{print NF}')
if [[ $ncol == 5 ]]; then
    zcat ${OUTPUT_RD_DIR}/${sample_name}.regions.bed.gz | cut -f1-3,5 | bgzip -c > ${OUTPUT_RD_DIR}/${sample_name}.cov.bed.gz
else
    echo "Error, please check the mosdepth output for (the number of columns) ${OUTPUT_RD_DIR}/${sample_name}.regions.bed.gz" 1>&2
    exit 1
fi

tabix -f -p bed ${OUTPUT_RD_DIR}/${sample_name}.cov.bed.gz

echo "Job ended on `hostname` at `date`"
