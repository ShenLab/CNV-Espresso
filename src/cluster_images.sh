#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=100:0:0
#$ -l h_vmem=18G
#S -l mem_free=18G
#$ -pe smp 1
set -oe pipefail

echo "Job started on `hostname` at `date`"
script_dir='/home/rt2776/cnv_espresso/src/'
rd_norm_dir=$1
ref_samples_dir=$2
cnv_file=$3
output_dir=$4

python ${script_dir}cnv_espresso.py images \
    --rd_norm_dir ${rd_norm_dir} \
    --ref_dir ${ref_samples_dir} \
    --cnv_list ${cnv_file} \
    --output ${output_dir} \
    --specific ${SGE_TASK_ID}

echo "Job ended on `hostname` at `date`"
