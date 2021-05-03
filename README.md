# _CNV Espresso_
#### A tool designed for validating **C**opy **N**umber **V**ariants from **E**xome **S**equencing **PRE**diction**S** in **S**ilic**O**

## How to run
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
    python /home/rt2776/cnv_espresso/src/cnv_espresso.py normalization \
        --windows ${project_dir}/windows.bed \ 
        --input /home/rt2776/1000GP/data/RD_clamms/NA12878.cov.bed.gz \ 
        --output ${project_dir}'/norm/'

    2. Multiple samples (input a simple list) 
    python /home/rt2776/cnv_espresso/src/cnv_espresso.py normalization \
        --windows ${project_dir}/windows.bed \ 
        --input_list ${project_dir}/sample_cov.list \
        --output ${project_dir}'/norm/'

    3. Multiple samples process by cluster
    windows_file='/home/rt2776/1000GP/cnv_espresso/windows_sort.bed'
    input_RD_list='/home/rt2776/1000GP/cnv_espresso/sample_rd.list'
    output_dir='/home/rt2776/1000GP/cnv_espresso/norm/'
    qsub -t 1-90 /home/rt2776/cnv_espresso/src/cluster_gc_norm.sh \
        ${windows_file} ${input_RD_list} ${output_dir} 

### Step 4. Calculate correlation matrix



### Step 5. Select reference samples
    ls /home/rt2776/1000GP/3_cnv_espresso_BI/norm/*.gz >/home/rt2776/1000GP/3_cnv_espresso_BI/sample_norm_rd.list

    '''
    /home/rt2776/cnv_espresso/src/select_reference.Rmd
    /home/rt2776/cnv_espresso/src/select_reference.ipynb (Preivous Prefered, now incorporated into cnv_espresso.py)
    '''

    python /home/rt2776/cnv_espresso/src/cnv_espresso.py reference \ 
        --project_path /home/rt2776/1000GP/3_cnv_espresso_BI/ \
        --norm_list /home/rt2776/1000GP/3_cnv_espresso_BI/sample_norm_rd.list \
        --num_ref 100 \
        --corr_threshold -1 

### Step 6. Generate images 

- Generate a single cnv
    RD_norm_dir='/home/rt2776/1000GP/3_cnv_espresso_BI/norm/'
    ref_samples_dir='/home/rt2776/1000GP/3_cnv_espresso_BI/ref_samples/'
    cnv_list='/home/rt2776/1000GP/3_cnv_espresso_BI/xhmm.xcnv' 
    output_dir='/home/rt2776/1000GP/3_cnv_espresso_BI/images_ref70_flanking/'
    corr_threshold=0.70
    flanking=True
    python /home/rt2776/cnv_espresso/src/generate_images.py \
        ${RD_norm_dir} ${ref_samples_dir} ${cnv_list} ${output_dir} ${corr_threshold} ${flanking} 1048 1037  1051

    qsub -t 1037-1051 /home/rt2776/cnv_espresso/src/cluster_generate_images.sh \
        ${RD_norm_dir} ${ref_samples_dir} ${cnv_list} ${output_dir} ${corr_threshold} ${flanking}

### images check and annotate file name into the cnv file for downstream analysis
    cnv_list='/home/rt2776/1000GP/3_cnv_espresso_BI/xhmm_NA12878.xcnv' 
    output_dir='/home/rt2776/1000GP/3_cnv_espresso_BI/images_ref70_flanking/'
    python /home/rt2776/cnv_espresso/src/generate_images_results_check_annotate.py \
        ${cnv_list} ${output_dir}


## predictions
To better analysis the prediction results, we need to annotate 1)batch info; 2)num_of_windows; 3)CNV frequency. 
Note: those CNV calls were generated from both WES and Array data. Therefore, it is normal for CNVs with 0 target.
This step is best run on a GPU machine. 

    cnv_w_image_file='/home/rt2776/1000GP/3_cnv_espresso_BI/images_ref70_flanking/xhmm_NA12878_withImagePath.csv'
    model_file='/home/rt2776/cnv_espresso/images_rare_3classes/data_backup/model_h5/rare_entire_cnv_MobileNet_v1_fine-tuning_3classes.h5'
    output_file='/home/rt2776/1000GP/3_cnv_espresso_BI/cnv_espresso_prediction_ref70_flanking.csv'
    python /home/rt2776/cnv_espresso/src/prediction.py \
        ${cnv_w_image_file} ${model_file} ${output_file}

