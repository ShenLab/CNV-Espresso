## Step 0. Configure the path

```bash
script_dir='/home/rt2776/cnv_espresso/src/'
project_dir='/home/rt2776/cnv_espresso/example/'
output_rd_dir=${project_dir}'/RD_data/'
target_file='/home/rt2776/source/capture_kits/20130108.exome.targets.bed'
reference_file='/home/rt2776/source/references/human_g1k_v37.fasta'
cd ${project_dir}
```

- SPARK WES2

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
- SPARK WES1 training model

  ```bash
  script_dir='/home/rt2776/cnv_espresso/src/'
  project_dir='/home/rt2776/cnv_espresso/project0_train_model/'
  output_rd_dir=${project_dir}'/RD_data/'
  target_file='/share/terra/SPARK_WES_2/SPARK_Freeze_Five_pVCF/resources/xgen_plus_spikein.b38.bed'
  reference_file='/share/terra/rsrc/hg38/ref/genome.hg38rg.fa'
  cd ${project_dir}
  
  # Step5
  RD_norm_dir=${project_dir}/norm_clamms/
  ref_samples_dir=${project_dir}/ref_samples/
  cnv_list=${project_dir}/spark_wes2_ultra_rare_cnvs.txt
  output_dir=${project_dir}

  # Step7
  model_file='/home/rt2776/cnv_espresso/model/MobileNet_v1_fine_tuning_3classes.h5'

  ```


### Step 2. Calculate Read depth for each sample

- Option 1. A single sample example

```bash
sample_name='NA12878'
bam_cram_file=/home/rt2776/1000GP/data/BAM/NA12878.mapped.ILLUMINA.bwa.CEU.exome.20121211.bam
```

