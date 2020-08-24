# CNV_Espresso
a tool designed for validating Copy Number Variants from Exome Sequencing PREdictionS in SilicO

## path
```
GATK_SW_DIR
TARGET_PROBES
REF_GENOME

SCRIPTS_DIR='/home/rt2776/CNV_Espresso/src/'
DATA_DIR='/home/rt2776/SPARK/CANOES/0-RC/spark_merged_RC/'
OUTPUT_DIR='/home/rt2776/CNV_Espresso/result/NormReadCountRatio/'

```

## Step 1. Calculate RD/RC for each sample in the given cohort

## Step 2. Extract GC content for each interval (Not essential)
```
java -Xmx2000m -Djava.io.tmpdir=${DATA_LOGS_DIR} \
                      -jar ${GATK_SW_DIR}GenomeAnalysisTK.jar \
                      -T GCContentByInterval \
                      -L ${TARGET_PROBES} \
                      -R ${REF_GENOME} \
                      -o ${DATA_DIR}gc.txt
```

## Step 3. Extract read count ratio for each sample
```
TEST_SAMPLE='SP0000203'
BATCH_NAME='spark1'
Rscript --vanilla ${SCRIPTS_DIR}extract_rc_ratio.R ${DATA_DIR} ${TEST_SAMPLE} ${BATCH_NAME} ${OUTPUT_DIR}   

bgzip ${OUTPUT_DIR}${TEST_SAMPLE}.txt
tabix -p bed ${OUTPUT_DIR}${TEST_SAMPLE}.txt.gz

```

