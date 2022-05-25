#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=100:0:0
#$ -l h_vmem=18G
#S -l mem_free=18G
#$ -pe smp 1

#SBATCH -o job.%j.out
#SBATCH -J image_human_view 
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 24G
#SBATCH --time=99:99:99

set -oe pipefail
echo "Job started on `hostname` at `date`"
#export PATH=$HOME/miniconda3/bin:$PATH #Note: add these if it only works in a specific env 
#source activate str

ID=${SGE_TASK_ID}${SLURM_ARRAY_TASK_ID}
cluster_input_file=$1

script_dir=$(head -n $ID $cluster_input_file | tail -n 1 |awk '{print $1}')
rd_norm_dir=$(head -n $ID $cluster_input_file | tail -n 1 |awk '{print $2}')
ref_samples_dir=$(head -n $ID $cluster_input_file | tail -n 1 |awk '{print $3}')
cnv_file=$(head -n $ID $cluster_input_file | tail -n 1 |awk '{print $4}')
vcf_file=$(head -n $ID $cluster_input_file | tail -n 1 |awk '{print $5}')
output_dir=$(head -n $ID $cluster_input_file | tail -n 1 |awk '{print $6}')
flanking=$(head -n $ID $cluster_input_file | tail -n 1 |awk '{print $7}')
suffix=$(head -n $ID $cluster_input_file | tail -n 1 |awk '{print $8}')
overwrite_img=$(head -n $ID $cluster_input_file | tail -n 1 |awk '{print $9}')

echo "ID:" $ID
echo "cluster_input_file:" ${cluster_input_file}
echo "script_dir:" ${script_dir}
echo "rd_norm_dir:" ${rd_norm_dir}
echo "ref_samples_dir:" ${ref_samples_dir}
echo "cnv_file:" ${cnv_file}
echo "vcf_file:" ${vcf_file}
echo "output_dir:" ${output_dir}
echo "flanking:" ${flanking}
echo "suffix:" ${suffix}
echo "overwrite_img:" ${overwrite_img}


python ${script_dir}cnv_espresso.py images_human_view \
    --rd_norm_dir ${rd_norm_dir} \
    --ref_dir ${ref_samples_dir} \
    --cnv_list ${cnv_file} \
    --vcf_file ${vcf_file} \
    --output ${output_dir} \
    --specific ${SGE_TASK_ID}${SLURM_ARRAY_TASK_ID} \
    --flanking ${flanking} \
    --suffix ${suffix} \
    --overwrite_img ${overwrite_img}

echo "Job ended on `hostname` at `date`"
