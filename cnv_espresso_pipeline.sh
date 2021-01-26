# _CNV Espresso_
#### A tool designed for validating **C**opy **N**umber **V**ariants from **E**xome **S**equencing **PRE**diction**S** in **S**ilic**O**

#### Step 0. Configure the path
```
GATK_SW_DIR
TARGET_PROBES
REF_GENOME

PROJECT_DIR='/home/rt2776/cnv_espresso/'
SCRIPTS_DIR='/home/rt2776/cnv_espresso/src/'
CNV_TOOLKIT_DIR='/home/rt2776/cnv_toolkit/scripts/'
DATA_DIR='/home/rt2776/SPARK/CANOES/0-RC/spark_merged_RC/'
OUTPUT_DIR='/home/rt2776/cnv_espresso/result/NormReadCountRatio/'

```

##--------------------------------------------------------------------------------------------------------------##
## Step 1. Calculate RC for each sample in the given cohort
# Approach 1: using Mosdepth
# Approach 2: by my script
#TODO

##--------------------------------------------------------------------------------------------------------------##
## Step 2. Calculate correlation matrix
#TODO
#- note: we need to test/determing the correlation matrix calculation should be set  before or after GC and z-score normalization.

##--------------------------------------------------------------------------------------------------------------##
## Step 3. Extract GC content for each interval (Not essential)
```
#TODO
java -Xmx2000m -Djava.io.tmpdir=${DATA_LOGS_DIR} \
                      -jar ${GATK_SW_DIR}GenomeAnalysisTK.jar \
                      -T GCContentByInterval \
                      -L ${TARGET_PROBES} \
                      -R ${REF_GENOME} \
                      -o ${DATA_DIR}gc.txt
```
##--------------------------------------------------------------------------------------------------------------##
## Step 4. GC normalization, fit the negative binomial model and estimate mu and sigma, zscore 
# just start with CLAMMS' windows.bed and Mosdepth cov.bed.gz files. We need to do it by ourselves later.
sort -k1,1V -k2,2n windows.bed >windows_sort.bed

python /home/rt2776/cnv_espresso/src/cnv_espresso.py normalization \
    --windows /home/rt2776/cnv_espresso/data/windows_sort.bed \
    --input /home/rt2776/cnv_espresso/data/ForCLAMMs/SP0023447.cov.bed.gz \
    --output /home/rt2776/cnv_espresso/data/norm

## alternatively, input a list with multiple samples 
python /home/rt2776/cnv_espresso/src/cnv_espresso.py normalization \
    --windows /home/rt2776/cnv_espresso/data/windows_sort.bed \
    --input_list /home/rt2776/cnv_espresso/data/ref/SP0000027.ref.list \
    --output /home/rt2776/cnv_espresso/data/norm

## Normalize all the samples of SPARK WES1
python /home/rt2776/cnv_espresso/src/cnv_espresso.py normalization \
    --windows /home/rt2776/cnv_espresso/data/windows_sort.bed \
    --input_list /home/rt2776/cnv_espresso/data/sample_rd.list \
    --output /home/rt2776/cnv_espresso/data/norm
### by cluster
windows_file='/home/rt2776/cnv_espresso/data/windows_sort.bed'
input_RD_list='/home/rt2776/cnv_espresso/data/sample_rd.list'
output_dir='/home/rt2776/cnv_espresso/data/norm'
qsub -tc 30 -t 1-28389 /home/rt2776/cnv_espresso/src/cluster_gc_norm.sh ${windows_file} ${input_RD_list} ${output_dir} 

### build index 
bgzip SP0023447.cov.bed.norm
tabix -p bed SP0023447.cov.bed.norm.gz

##--------------------------------------------------------------------------------------------------------------##
## Step 5. Select reference samples
sample_list='/home/rt2776/cnv_espresso/source/spark_sample_27270.txt'
sample_rd_dir='/home/rt2776/cnv_espresso/data/norm'

### Open with Rstudio
/home/rt2776/cnv_espresso/src/select_reference.Rmd

##--------------------------------------------------------------------------------------------------------------##
## Step 6. Generate fixed-size windows for the CNV predicted regions
### the training cnvs should be rare cnvs, otherwise it will confuse the model
python ~/cnv_toolkit/scripts/5_annotation.py cnvfrequency_given_cohort \
    --input /home/rt2776/cnv_espresso/training_set/training_set_false.txt \
    --given_cohort /home/rt2776/cnv_control/Control_SPARK27_unaffected_parents.txt \
    --output /home/rt2776/cnv_espresso/training_set/training_set_false_af.txt 

data_path='/home/rt2776/cnv_espresso/training_set/'
qsub -t 1-589 /home/rt2776/cnv_toolkit/scripts/cluster_af_give_cohort.sh \
     ${data_path}training_set_true.txt \
     /home/rt2776/cnv_control/Control_SPARK27_unaffected_parents.txt \
     ${data_path}af_true/training_set_true_af.txt

cat ${data_path}af_true/*core* >${data_path}training_set_true_af.txt
### Filtering af >1% in cohort
NUM_SAMPLES=`cat /home/rt2776/cnv_control/Control_SPARK27_unaffected_parents.txt |awk '{print $5}'|sort|uniq|wc -l`
THRESH_NUM=`echo "0.01 * $NUM_SAMPLES"|bc|awk '{x = $1; if (x != int(x)) {x = int(x)+1} print (x-1)}'`
cat /home/rt2776/cnv_espresso/training_set/tmp/training_set_true_af.txt |awk -F '\t' \
    '{if($19 == "Num_Carriers(inGivenCohort)" || $19 <= '$THRESH_NUM') print $0}' >${data_path}training_set_true_rare.txt

cat /home/rt2776/cnv_espresso/training_set/tmp/training_set_false_af.txt |awk -F '\t' \
    '{if($19 == "Num_Carriers(inGivenCohort)" || $19 <= '$THRESH_NUM') print $0}' >${data_path}training_set_false_rare.txt

### Focus on 27k samples
cat training_set_true_rare.txt |grep -w -f ../source/spark_sample_27270.txt > training_set_true_rare_27k.txt
cat training_set_false_rare.txt |grep -w -f ../source/spark_sample_27270.txt > training_set_false_rare_27k.txt

### debuging by jupyter-lab
/home/rt2776/cnv_espresso/src/generate_images.ipynb

### illustrate single cnv
RD_norm_dir='/home/rt2776/cnv_espresso/data/norm/'
ref_samples_dir='/home/rt2776/cnv_espresso/reference_samples/'
python /home/rt2776/cnv_espresso/src/generate_images.py \
    ${RD_norm_dir} ${ref_samples_dir} \
    /home/rt2776/cnv_espresso/training_set/training_set_true_rare_27k.txt \
    /home/rt2776/cnv_espresso/images_rare/ 2288


RD_norm_dir='/home/rt2776/cnv_espresso/data/norm/'
ref_samples_dir='/home/rt2776/cnv_espresso/reference_samples/'
qsub -t 1-6927 /home/rt2776/cnv_espresso/src/cluster_generate_images.sh \
    ${RD_norm_dir} ${ref_samples_dir} \
    /home/rt2776/cnv_espresso/training_set/training_set_false_rare.txt \
    /home/rt2776/cnv_espresso/images_rare/ 

qsub -t 1-16371 /home/rt2776/cnv_espresso/src/cluster_generate_images.sh \
    ${RD_norm_dir} ${ref_samples_dir} \
    /home/rt2776/cnv_espresso/training_set/training_set_true_rare.txt \
    /home/rt2776/cnv_espresso/images_rare/ 

### output images check. check whether the number of images equal to the number of CNVs in given.
### conclustion: the missing CNVs should not belong to SPARK WES1 27270 samples.
RD_norm_dir='/home/rt2776/cnv_espresso/data/norm/'
ref_samples_dir='/home/rt2776/cnv_espresso/reference_samples/'
python /home/rt2776/cnv_espresso/src/generate_images_results_check.py \
    ${RD_norm_dir} ${ref_samples_dir} \
    /home/rt2776/cnv_espresso/training_set/training_set_true_rare_27k.txt \
    /home/rt2776/cnv_espresso/images_rare/
   
##--------------------------------------------------------------------------------------------------------------##
## Step 7. Deep Learning model



##--------------------------------------------------------------------------------------------------------------##
# Previous scripts in below

#### Step 3a. Extract read count ratio for each sample
```
TEST_SAMPLE='SP0000203'
BATCH_NAME='spark1'
Rscript --vanilla ${SCRIPTS_DIR}extract_rc_ratio.R ${DATA_DIR} ${TEST_SAMPLE} ${BATCH_NAME} ${OUTPUT_DIR}   

bgzip ${OUTPUT_DIR}${TEST_SAMPLE}.txt
tabix -p bed ${OUTPUT_DIR}${TEST_SAMPLE}.txt.gz

```
#### Step 3a. [Cluster] Extract read count ratio for each sample by cluster
SAMPLE_W_BATCH_LIST=/home/rt2776/SPARK/CANOES/0-RC/spark_27270_by10Groups/spark_w_batch.list
#try to avoid "OMP: Error #34: System unable to allocate necessary resources for OMP thread:"
export OMP_NUM_THREADS=1
SAMPLE_W_BATCH_LIST=/home/rt2776/SPARK/CANOES/0-RC/spark_27270_by10Groups/spark_w_batch_offspring.list
qsub -tc 36 -t 1-12744 ${SCRIPTS_DIR}extract_rc_ratio_cluster.sh ${SCRIPTS_DIR} ${DATA_DIR} ${SAMPLE_W_BATCH_LIST} ${OUTPUT_DIR}

#### Step 3b. Checking results to make sure every sample has been extracted.
sh ${SCRIPTS_DIR}extract_rc_ratio_sample_check.sh ${SAMPLE_W_BATCH_LIST} ${OUTPUT_DIR}

#### Step 4. Extract Read Count info for training set
python ${CNV_TOOLKIT_DIR}5_annotation.py rc_ratio \
    --input /home/rt2776/CNV_Espresso/result/training_set_final_all_clean_af.cnv.bz2 \
    --rc_ratio /home/rt2776/CNV_Espresso/result/NormReadCountRatio \
    --output /home/rt2776/CNV_Espresso/result/training_set_w_rc.cnv

#### Step 5. Split Read Count info in a CNV into multiple same-size windows
python ${SCRIPTS_DIR}cnv_espresso.py split \
    --input ${PROJECT_DIR}result/training_set_w_rc_clean.cnv \
    --output ${PROJECT_DIR}result/training_set_w_rc_equal_win.cnv \
    --num_target 3
