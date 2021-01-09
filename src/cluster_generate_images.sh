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
cnv_file=$1
python /home/rt2776/cnv_espresso/src/generate_images.py ${cnv_file} ${SGE_TASK_ID}

echo "Job ended on `hostname` at `date`"
