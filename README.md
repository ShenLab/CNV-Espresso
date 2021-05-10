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
    cnv_list=${project_dir}/xhmm.xcnv
    output_dir=${project_dir}

    1. Generate images
    python ${script_dir}cnv_espresso.py images \
        --rd_norm_dir ${RD_norm_dir} \
        --ref_dir ${ref_samples_dir} \
        --cnv_list ${cnv_file} \
        --output ${output_dir} \
        --specific 10 

    2. Generate images by cluster 
    ^ Note: you need to modify the script path in the `cluster_generate_images.sh` at first

    qsub -t 1-1505 ${script_dir}cluster_images.sh \
        ${RD_norm_dir} ${ref_samples_dir} ${cnv_list} ${output_dir} 


### Step 6. Training 
    In general, you can directly use our pretrained CNN model(/URL) for your in silico validation (Skip step 6). However, if you have a bunch of validated or confirmed CNVs, you can also train the CNN model from scratch. If so, please follow the tutorial below:

    1. Please use `images` function in `cnv_espresso.py` as *step5* to generate images for your prepared true deletion, true duplication, false deletion and false duplication. Note that false deletion and false duplication will be treated as diploid together in the downstream steps.

    2. train the model by the following command:
    true_del_img=${project_dir}/train/true_del.list
    true_dup_img=${project_dir}/train/true_dup.list
    false_del_img=${project_dir}/train/false_del.list
    false_dup_img=${project_dir}/train/false_dup.list

    python ${script_dir}cnv_espresso.py train \
        --true_del ${true_del_img} \
        --true_dup ${true_dup_img} \
        --false_del ${false_del_img} \
        --false_dup ${false_dup_img} \
        --use_gpu True \
        --output ${output_dir}
    
    Alternatively, we also prepared a jupyter notebook(/URL) for tracking and debugging the entire training process.

### Step 7. Validating CNV predictions in silico 
    cnv_w_img_file=${project_dir}/cnv_info_w_img.csv
    model_file='/home/rt2776/cnv_espresso/model/rare_entire_cnv_MobileNet_v1_fine-tuning_3classes.h5'
    output_file=${project_dir}/cnv_espresso_prediction_cpu.csv

    python ${script_dir}cnv_espresso.py predict \
        --cnv_list ${cnv_w_img_file} \
        --model ${model_file} \
        --output ${output_file} \
        --use_gpu False

## [TODO] Other functions
    1. plot read depth before and after GC normlization

    2. Merge results from multiple CNV callers

