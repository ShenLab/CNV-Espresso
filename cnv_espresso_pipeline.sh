# _CNV Espresso_
#### A tool designed for validating **C**opy **N**umber **V**ariants from **E**xome **S**equencing **PRE**diction**S** in **S**ilic**O**

#### Step 0. Configure the path
```
GATK_SW_DIR
TARGET_PROBES
REF_GENOME

SCRIPTS_DIR='/home/rt2776/CNV_Espresso/src/'
DATA_DIR='/home/rt2776/SPARK/CANOES/0-RC/spark_merged_RC/'
OUTPUT_DIR='/home/rt2776/CNV_Espresso/result/NormReadCountRatio/'

```

#### Step 1. Calculate RD/RC for each sample in the given cohort

#### Step 2. Extract GC content for each interval (Not essential)
```
java -Xmx2000m -Djava.io.tmpdir=${DATA_LOGS_DIR} \
                      -jar ${GATK_SW_DIR}GenomeAnalysisTK.jar \
                      -T GCContentByInterval \
                      -L ${TARGET_PROBES} \
                      -R ${REF_GENOME} \
                      -o ${DATA_DIR}gc.txt
```

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


#### Step 5. Extract Read Count info for testing set


#### Step 6. Build the Deep Learning model 
