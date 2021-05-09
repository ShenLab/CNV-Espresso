# _CNV Espresso_
#### A tool designed for validating **C**opy **N**umber **V**ariants from **E**xome **S**equencing **PRE**diction**S** in **S**ilic**O**

## Usage
### Step 0. Configure the path
    script_dir='/home/rt2776/cnv_espresso/src/'
    project_dir='/home/rt2776/cnv_espresso/project4_method_development'
    target_file='/home/rt2776/source/capture_kits/20130108.exome.targets.bed'
    reference_file='/home/rt2776/source/references/human_g1k_v37.fasta'
    cd ${project_dir}

### Step 1. Calculate Read depth for each sample
    Mosdepth

### Step 2. Split large capture region into small windows and annotate GC content 
    python ${script_dir}cnv_espresso.py windows \
        --target ${target_file} \
        --ref    ${reference_file} \
        --output ${project_dir}

### Step 3. GC normalization
    1. Single sample
    python ${script_dir}cnv_espresso.py normalization \
        --windows ${project_dir}/windows.bed \ 
        --input /home/rt2776/1000GP/data/RD_clamms/NA12878.cov.bed.gz \ 
        --output ${project_dir}'/norm/'

    2. Multiple samples (input a simple list) 
    python ${script_dir}cnv_espresso.py normalization \
        --windows ${project_dir}/windows.bed \ 
        --input_list ${project_dir}/sample_raw_rd.list \
        --output ${project_dir}'/norm/'

    3. Multiple samples process by cluster 
    ^ TODO: need to test.
    windows_file='/home/rt2776/1000GP/cnv_espresso/windows_sort.bed'
    input_RD_list='/home/rt2776/1000GP/cnv_espresso/sample_raw_rd.list'
    output_dir='/home/rt2776/1000GP/cnv_espresso/norm/'
    qsub -t 1-90 /home/rt2776/cnv_espresso/src/cluster_gc_norm.sh \
        ${windows_file} ${input_RD_list} ${output_dir} 

### Step 4. Select reference samples
    ls ${project_dir}/norm/*.gz > ${project_dir}/sample_norm_rd.list

    python ${script_dir}cnv_espresso.py reference \ 
        --project_dir ${project_dir} \
        --norm_list ${project_dir}/sample_norm_rd.list \
        --num_ref 100 \
        --corr_threshold -1 

    ^ TODO: This step can show time. it is helpful for other functions.

### Step 5. Generate images 
    RD_norm_dir=${project_dir}/norm/
    ref_samples_dir=${project_dir}/ref_samples/
    cnv_file=${project_dir}/xhmm.xcnv
    output_dir=${project_dir}

    1. Generate images
    python ${script_dir}cnv_espresso.py images \
        --rd_norm_dir ${RD_norm_dir} \
        --ref_dir ${ref_samples_dir} \
        --cnv_list ${cnv_file} \
        --output ${output_dir} \
        --specific 1048 

    2. Generate images by cluster 
    ^ Note: you need to modify the script path in the `cluster_generate_images.sh` at first

    qsub -t 1-1505 ${script_dir}cluster_images.sh \
        ${RD_norm_dir} ${ref_samples_dir} ${cnv_list} ${output_dir} 

    qsub -t 20-35 ${script_dir}cluster_images.sh \
        ${RD_norm_dir} ${ref_samples_dir} ${cnv_list} ${output_dir} 

### Step 6. Training 

### Step 7. Validating CNV predictions in silico 
    cnv_w_img_file=${project_dir}/cnv_info_w_img.csv
    model_file='/home/rt2776/cnv_espresso/model/rare_entire_cnv_MobileNet_v1_fine-tuning_3classes.h5'
    output_file=${project_dir}/cnv_espresso_prediction.csv

    python ${script_dir}cnv_espresso.py predict \
        --cnv_list ${cnv_w_img_file} \
        --model ${model_file} \
        --output ${output_file} \
        --use_gpu False

## [TBD] Other functions
### plot read depth before and after GC normlization

### Merge results from multiple CNV callers

