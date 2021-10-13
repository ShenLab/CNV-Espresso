## SPARK WES2

  ```bash
  script_dir='/home/rt2776/cnv_espresso/src/'
  project_dir='/home/rt2776/cnv_espresso/project5_spark_wes2/'
  output_rd_dir='/home/rt2776/spark_wes2/data/RD_clamms_16004/'
  target_file='/share/terra/SPARK_WES_2/SPARK_Freeze_Five_pVCF/resources/xgen_plus_spikein.b38.bed'
  reference_file='/share/data/resources/hg38/references/genome.hg38rg.fa.gz'
  cd ${project_dir}
  
  # Step5
  RD_norm_dir=${project_dir}/norm_clamms/
  ref_samples_dir=${project_dir}/ref_samples/
  cnv_list=${project_dir}/spark_wes2_ultra_rare_cnvs.txt
  output_dir=${project_dir}

  # Step7
  model_file='/home/rt2776/cnv_espresso/model/MobileNet_v1_fine_tuning_3classes.h5'
  ```

## SPARK WES1 training model

  ```bash
  script_dir='/home/rt2776/cnv_espresso/src/'
  project_dir='/home/rt2776/cnv_espresso/project0_train_model/'
  output_rd_dir=${project_dir}'/RD_data/'
  target_file='/share/terra/SPARK_WES_2/SPARK_Freeze_Five_pVCF/resources/xgen_plus_spikein.b38.bed'
  reference_file='/share/terra/rsrc/hg38/ref/genome.hg38rg.fa'
  cd ${project_dir}
  
  # Step 2. Calculate Read depth for each sample
  ## fix a sample 
    sample_name='SP0000140'
    bam_cram_file='/share/terra/xz2680/SPARK_30K/wes.DeDupedCRAMs/SP0000140.cram'
    windows_file=${project_dir}/windows.bed
    mkdir -p ${output_rd_dir}

    sample_name='SP0065987'
    bam_cram_file='/share/terra/xz2680/SPARK_30K/wes.DeDupedCRAMs/SP0065987.cram'
    windows_file=${project_dir}/windows.bed
    mkdir -p ${output_rd_dir}

    mosdepth -n --fasta ${reference_file} --by ${windows_file} --mapq 30 ${output_rd_dir}/${sample_name} ${bam_cram_file}
    zcat ${output_rd_dir}/${sample_name}.regions.bed.gz | cut -f1-3,5 | bgzip -c > ${output_rd_dir}/${sample_name}.cov.bed.gz
    tabix -f -p bed ${output_rd_dir}/${sample_name}.cov.bed.gz

    ### Step 3. GC normalization
    python ${script_dir}cnv_espresso.py normalization \
        --windows ${project_dir}/windows.bed \
        --input   ${output_rd_dir}/SP0000140.cov.bed.gz \
        --output  ${project_dir}/norm/

  # Step5
  RD_norm_dir=${project_dir}/norm/
  ref_samples_dir=${project_dir}/ref_samples/

  cnv_list=${project_dir}/images_rare_3classes/0-entire_cnv_file_list/false_del_image_info.list
  output_dir=${project_dir}/logDiffCumX_logY_false_del/

  cnv_list=${project_dir}/images_rare_3classes/0-entire_cnv_file_list/false_dup_image_info.list
  output_dir=${project_dir}/logDiffCumX_logY_false_dup/

  cnv_list=${project_dir}/images_rare_3classes/0-entire_cnv_file_list/true_del_image_info.list
  output_dir=${project_dir}/logDiffCumX_logY_true_del/

  cnv_list=${project_dir}/images_rare_3classes/0-entire_cnv_file_list/true_dup_image_info.list
  output_dir=${project_dir}/logDiffCumX_logY_true_dup/

    python ${script_dir}cnv_espresso.py images \
        --rd_norm_dir ${RD_norm_dir} \
        --ref_dir     ${ref_samples_dir} \
        --cnv_list    ${cnv_list} \
        --output      ${output_dir}


  # Step7
  model_file='/home/rt2776/cnv_espresso/model/MobileNet_v1_fine_tuning_3classes.h5'

  ```

## project7 spark wes3 wes+wgs
  ```bash
  ## For slurm
  script_dir='/home/nas-0-1/nova.home/rt2776/cnv_espresso/src/'
  project_dir='/home/nas-0-1/nova.home/rt2776/cnv_espresso/project7_spark_wes3_wes_wgs/02b_clamms_cnvEspresso/'

  ## For sge
  script_dir='/home/rt2776/cnv_espresso/src/'
  project_dir='/home/rt2776/cnv_espresso/project7_spark_wes3_wes_wgs/02b_clamms_cnvEspresso/'

  output_rd_dir=${project_dir}'/RD_data/'
  target_file='/share/terra/SPARK_WES_2/SPARK_Freeze_Five_pVCF/resources/xgen_plus_spikein.b38.bed'
  reference_file='/share/terra/rsrc/hg38/ref/genome.hg38rg.fa'
  mkdir -p ${project_dir}
  cd ${project_dir}
  
  # Run the pipeline in `Readme.md` for Step1 to Step7.

  # Step5
  cnv_list=${project_dir}/clamms_high_confidence.cnv.gsd_anno_af.txt

  # Step7
  model_file='/home/nas-0-1/nova.home/rt2776/cnv_espresso/model/MobileNet_v1_fine-tuning_3classes_log_transformed.h5'

  ```
