# _CNV-Espresso_
#### A tool designed for validating **C**opy **N**umber **V**ariants from **E**xome **S**equencing **PRE**diction**S** in **S**ilic**O**


## Dependencies

Python with the following required libraries: tensorflow, keras, sklearn, numpy, pandas, matplotlib, seaborn, Pillow, pysam and lockfile

Mosdepth: https://github.com/brentp/mosdepth

Tabix: http://www.htslib.org/doc/tabix.html (optional)

## Usage

### Step 0. Configure the path
```bash
script_dir='/path/to/cnv_espresso/src/'
project_dir='/path/to/cnv_espresso/example/'
output_rd_dir=${project_dir}'/RD_data/'
target_file='/path/to/exome.targets.bed'
reference_file='/path/to/reference.fasta'
cd ${project_dir}
```

### Step 1. Split large capture region into small windows and annotate with GC content 
```bash
python ${script_dir}cnv_espresso.py windows \
    --target ${target_file} \
    --ref    ${reference_file} \
    --output ${project_dir}
```

### Step 2. Calculate Read depth for each sample

- Option 1. A single sample example

```bash
sample_name='example_sample'
bam_cram_file=/path/to/example_sample.cram
windows_file=${project_dir}/windows.bed
mkdir -p ${output_rd_dir}

mosdepth -n --fasta ${reference_file} --by ${windows_file} --mapq 30 ${output_rd_dir}/${sample_name} ${bam_cram_file}
zcat ${output_rd_dir}/${sample_name}.regions.bed.gz | cut -f1-3,5 | bgzip -c > ${output_rd_dir}/${sample_name}.cov.bed.gz
tabix -f -p bed ${output_rd_dir}/${sample_name}.cov.bed.gz
```

This is a relatively time-consuming step, however, you can use the following command to accelerate the speed if you have a cluster.

- Option 2. Multiple samples (by cluster)

```bash
bam_cram_file_path_list=${project_dir}/sample_rd_file_list.txt
windows_file=${project_dir}/windows.bed

# Suppose we want to calculate RD for 1000 samples
qsub -t 1-1000 ${script_dir}cluster_mosdepth.sh \
    ${bam_cram_file_path_list} ${reference_file} ${windows_file} ${output_rd_dir}      
```

Collect all the coverage files for downstream usage. 

```bash
ls ${output_rd_dir}*.cov.bed.gz > ${project_dir}/sample_raw_rd.txt 
```

### Step 3. GC normalization

- Option 1. Single sample

```bash
python ${script_dir}cnv_espresso.py normalization \
    --windows ${project_dir}/windows.bed \
    --input   ${output_rd_dir}/example_sample.cov.bed.gz \
    --output  ${project_dir}/norm/
```

- Option 2. Multiple samples (input a simple list) 

```bash
python ${script_dir}cnv_espresso.py normalization \
    --windows    ${project_dir}/windows.bed \
    --input_list ${project_dir}/sample_raw_rd.txt \
    --output     ${project_dir}/norm/
```

- Option 3. Multiple samples (by cluster) 

```bash
windows_file=${project_dir}/windows.bed
RD_list=${project_dir}/sample_raw_rd.list
output_dir=${project_dir}/norm/

# Suppose we want to process 1000 samples
qsub -t 1-1000 ${script_dir}cluster_gc_norm.sh \
    ${windows_file} ${RD_list} ${output_dir}
```

### Step 4. Select reference samples

```bash
ls ${project_dir}/norm/*.gz > ${project_dir}/sample_norm_rd.list
```

```bash
python ${script_dir}cnv_espresso.py reference \
    --project_dir ${project_dir} \
    --norm_list   ${project_dir}/sample_norm_rd.list \
    --num_ref     100 \
    --corr_threshold -1 
```

### Step 5. Generate images 

Here, we will take other CNV caller's output (e.g. `xhmm.xcnv`) as our input and we will use the following command to encode those CNV predictions into images.

```bash
RD_norm_dir=${project_dir}/norm/
ref_samples_dir=${project_dir}/ref_samples/
cnv_list=${project_dir}/xhmm.xcnv # other CNV caller's output
output_dir=${project_dir}
```

- Option 1. Generate images by single thread

```bash
python ${script_dir}cnv_espresso.py images \
    --rd_norm_dir ${RD_norm_dir} \
    --ref_dir     ${ref_samples_dir} \
    --cnv_list    ${cnv_list} \
    --output      ${output_dir} 
```

- Option 2. Generate images by cluster
  Note: please modify the script path in the `cluster_images.sh` at first.

```bash
qsub -t 1-1000 ${script_dir}cluster_images.sh \
    ${RD_norm_dir} ${ref_samples_dir} ${cnv_list} ${output_dir} 
```

A few example images are located [here](https://github.com/ShenLab/CNV-Espresso/tree/main/example/images). 

### Step 6. Training 

In general, you can directly use our pretrained CNN model ([MobileNet_v1_fine_tuning_3classes.h5](https://github.com/ShenLab/CNV-Espresso/blob/main/model/MobileNet_v1_fine_tuning_3classes.h5)) for your *in silico* validation (Skip step 6). However, if you have a bunch of validated or confirmed CNVs, you can also train the CNN model from scratch. If so, please follow the tutorial below:

1. Please use `images` function in `cnv_espresso.py` as **Step5** to generate images for your prepared true deletion, true duplication, false deletion and false duplication. Note that false deletion and false duplication will be treated as diploid together in the downstream steps.

2. We recommend using a GPU server to handle this step. You can train your specific model by the following command:

```bash
true_del_img=${project_dir}/train/true_del.list
true_dup_img=${project_dir}/train/true_dup.list
false_del_img=${project_dir}/train/false_del.list
false_dup_img=${project_dir}/train/false_dup.list
output_dir=${project_dir}/model/

python ${script_dir}cnv_espresso.py train \
    --true_del  ${true_del_img} \
    --true_dup  ${true_dup_img} \
    --false_del ${false_del_img} \
    --false_dup ${false_dup_img} \
    --use_gpu   True \
    --output    ${output_dir}
```

Alternatively, we also prepared a jupyter notebook (**[train.ipynb](https://github.com/ShenLab/CNV-Espresso/blob/main/src/train.ipynb)**) for tracking and debugging the entire training process.

### Step 7. Validating CNV predictions *in silico* 

```bash
cnv_w_img_file=${project_dir}/cnv_info_w_img.csv
model_file='/path/to/model/MobileNet_v1_fine_tuning_3classes.h5'
output_file=${project_dir}/cnv_espresso_prediction.csv

python ${script_dir}cnv_espresso.py predict \
    --cnv_list ${cnv_w_img_file} \
    --model    ${model_file} \
    --output   ${output_file} \
    --use_gpu  False
```

## Utilities
We will release the following auxiliary functions in the near future.
- Plot read depth signal before and after GC normalization
- Merge results from multiple CNV callers

## Contact

Renjie Tan (rt2776 at cumc.columbia.edu) or Yufeng Shen (ys2411 at cumc.columbia.edu)

#### [Shen Lab](http://www.columbia.edu/~ys2411/)

