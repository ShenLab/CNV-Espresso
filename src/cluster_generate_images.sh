#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=100:0:0
#$ -l h_vmem=18G
#S -l mem_free=18G
#$ -pe smp 1
#$ -o '/home/rt2776/cnv_espresso/oe/'
#$ -e '/home/rt2776/cnv_espresso/oe/'

echo "Job started on `hostname` at `date`"
RD_norm_dir=$1
ref_samples_dir=$2
cnv_file=$3
output_dir=$4
corr_threshold=$5

python /home/rt2776/cnv_espresso/src/generate_images.py \
    ${RD_norm_dir} \
    ${ref_samples_dir} \
    ${cnv_file} \
    ${output_dir} \
    ${corr_threshold} \
    ${SGE_TASK_ID}

echo "Job ended on `hostname` at `date`"
