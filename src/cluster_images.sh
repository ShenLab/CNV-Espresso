#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=100:0:0
#$ -l h_vmem=18G
#S -l mem_free=18G
#$ -pe smp 1

#SBATCH -o ./logs/job.%j.out
#SBATCH -J image 
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 18G
#SBATCH --time=99:99:99

set -oe pipefail
echo "Job started on `hostname` at `date`"
script_dir=$1
rd_norm_dir=$2
ref_samples_dir=$3
cnv_file=$4
output_dir=$5
overwrite_img=$6

python ${script_dir}cnv_espresso.py images \
    --rd_norm_dir ${rd_norm_dir} \
    --ref_dir ${ref_samples_dir} \
    --cnv_list ${cnv_file} \
    --output ${output_dir} \
    --specific ${SGE_TASK_ID}${SLURM_ARRAY_TASK_ID} \
    --overwrite_img ${overwrite_img}

echo "Job ended on `hostname` at `date`"
