########################################################################
# Script : extract_bp_coverage.sh                                      #
# Author : Renjie Tan                                                  #
# Date   : 5/12/2019                                                    #
#                                                                      #
# This script calculates basepair level coverage information for       #
# each sample. This information is used later to resolve breakpoints   #
# of concordant CNVs.                                                  #
#                                                                      #
########################################################################
#!/bin/bash

echo "Job started on `hostname` at `date`"
sample_file='/home/rt2776/CN_Learn/source/spark_cram_list_20190729_merged.txt'
DATA_BPCOV_DIR='/share/terra/rt2776/SPARK/CN_Learn/bp_coverage_dir/'

num=0
IFS=$'\n' #use IFS variable to specific you want a newline as the field separator
for sample_reader in $(tac ${sample_file});
do
    dt=$(date '+%m/%d/%Y %H:%M:%S')
    #let "num++"
    #echo "[$dt] ${num} Processing ${sample_reader} .............."
    cram_file=$(echo ${sample_reader} | cut -d " " -f 1)
    sample_name=$(echo ${sample_reader} | cut -d " " -f 2)
    if [ ! -f ${DATA_BPCOV_DIR}${sample_name}.bpcov.bed.gz ];
    then
        let "num++"
        echo "[$dt] ${num} Extract bp coverage for sample:${sample_name}" 
        /home/rt2776/softwares/bedtools2/bin/genomeCoverageBed -ibam ${cram_file} -bga \
                             > ${DATA_BPCOV_DIR}${sample_name}.bpcov.bed

        bgzip ${DATA_BPCOV_DIR}${sample_name}.bpcov.bed
    fi
done
echo "Job ended on `hostname` at `date`"
