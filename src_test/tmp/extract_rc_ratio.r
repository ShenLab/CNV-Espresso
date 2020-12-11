###################################################################################
# Author: Renjie Tan                                                          
# Date  : 8/24/2020                                        
# This script is supplied as part of the CNV-Espresso pipeline. 
# It has been modified from CANOES pipeline to calculate the read count ratio.
# The key algorithm includes:
#   1) Normalize each sample's total read count is equal.
#   2) Select reference samples by covariantion
#   3) Weight the reference samples by non-negative lease square.
#   4) Calculate the reference read count by the weighted reference samples
###################################################################################

############################################################################
# STEP 1: Fetch the directory locations and file names from the shell script
############################################################################
args <- commandArgs(trailingOnly=TRUE)
data_path    <- args[1]
case_sample  <- args[2] # eg. "SP0000203"
batch_name   <- args[3]
output_path  <- args[4]

# sample list
sample_name_file <- paste(data_path,'merged_sample_list_', batch_name,'.list',sep='')

# output file
output_path <- paste(output_path, batch_name,'/',sep='')
output_file_name <- paste(output_path, 'NormRCRatio_',test_sample.name,'.txt',sep='')

######################################################
# STEP 2: Read the read count file and GC content file
######################################################
# gc file
gc_file <- paste(data_path, 'gc.txt.bz2',sep='')
file.exists(gc_file)
gc <- read.table(bzfile(gc_file))$V2

# input file
input_file_name <- paste(data_path, 'cons_canoes_reads_',batch_name,'.bz2',sep='') 
input_file_RDate <- paste(input_file_name,'.RData',sep='') 
if(file.exists(input_file_RDate)){
  load(input_file_RDate)
}else{
  file.exists(input_file_name)
  canoes.reads <- read.table(bzfile(input_file_name))
  save(canoes.reads,file = input_file_RDate)
}
canoes.reads <- canoes.reads[which(canoes.reads$target!='NA'),]

########################################################
# STEP 3: Assign column names and merge columns together
########################################################
num_cols <- dim(canoes.reads)[2]
num_samples <- num_cols - 3
sample.names <- read.table(file=sample_name_file)
sample.names <- sub(".tar.bz2","",sample.names[,1])
sample.names <- sub(".cram","",sample.names)
sample.names <- sub(".bam","",sample.names)
names(canoes.reads) <- c("chromosome", "start", "end", sample.names)
target <- seq(1, nrow(canoes.reads))
canoes.reads <- cbind(target, canoes.reads[,1:3], gc, canoes.reads[,-c(1,2,3)])

######################################################################################################
# STEP 4: Remove sex chromosomes and Force an object to belong to a class
######################################################################################################
if(class(canoes.reads$chromosome) != "numeric"){
  canoes.reads$chromosome <- as.character(canoes.reads$chromosome) #RT modifided
}

if (length(setdiff(names(canoes.reads)[1:5], c("target", "chromosome", "start", "end", "gc"))) > 0){
  stop("First five columns of counts matrix must be target, chromosome, start, end, gc")
}
if (length(setdiff(unique(canoes.reads$chromosome), seq(1:22))) > 0) {
  # remove sex chromosomes
  cat("Trying to remove sex chromosomes and 'chr' prefixes\n")
  canoes.reads <- subset(canoes.reads, !chromosome %in% c("chrX", "chrY", "X", "Y"))
  if (sum(grepl("chr", canoes.reads$chromosome))==length(canoes.reads$chromosome)){
    canoes.reads$chromosome <- gsub("chr", "", canoes.reads$chromosome)
  }
  canoes.reads$chromosome <- as.numeric(canoes.reads$chromosome)
  if (length(setdiff(unique(canoes.reads$chromosome), seq(1:22))) > 0) 
    stop("chromosome must take value in range 1-22 (support for sex chromosomes to come)")
}

if(class(canoes.reads$chromosome) != "numeric"){
  canoes.reads$chromosome <- as.numeric(as.character(canoes.reads$chromosome)) #RT modifided
}

######################################################################################################
# STEP 5: calculate covariance of read count across samples
######################################################################################################
file_name <- strsplit(basename(input_file_name),"\\.")[[1]][1]
cov_file  <- paste(data_path,'cov.',file_name,sep='')
if (file.exists(cov_file)){
  load(file=cov_file)
}
if (!exists('covariance')){
  cat("Please calculate covariance of read count across samples at first.")
  cat(paste('[',format(Sys.time()),']'," calculate covariance of read count across samples\n",sep=''))#
  system.time(covariance <- cor(canoes.reads[, sample.names], canoes.reads[, sample.names]))
  save(covariance, file=paste(data_path, 'cov.', file_name, sep=''))
}

######################################################
# STEP 6: Load R script and process each column/sample
######################################################
if (!case_sample %in% names(canoes.reads)){stop("No column for sample ", case_sample, " in counts matrix")}
if (length(setdiff(names(canoes.reads)[1:5], c("target", "chromosome", "start", "end", "gc"))) > 0){
  stop("First five columns of counts matrix must be target, chromosome, start, end, gc")
}

if(class(canoes.reads$chromosome) != "numeric"){
  canoes.reads$chromosome <- as.numeric(as.character(canoes.reads$chromosome)) #RT modifided
}

library(plyr)
canoes.reads <- arrange(canoes.reads, chromosome, start)
case_samples <- colnames(canoes.reads)[-seq(1,5)]

# find mean coverage of probes
cat(paste('[',format(Sys.time()),']'," find mean coverage of probes\n",sep=''))
mean.counts <- mean(apply(counts[, case_samples], 2, mean))

# normalize counts; round so we can use negative binomial
cat(paste('[',format(Sys.time()),']'," normalize counts; round so we can use negative binomial\n",sep=''))
counts[, case_samples] <- apply(counts[, case_samples], 2, 
                                function(x, mean.counts) 
                                  round(x * mean.counts / mean(x)), mean.counts)

# calculate covariance of read count across samples
reference.samples <- setdiff(case_samples, case_sample)
covariances <- covariance[case_sample, reference.samples]
reference.samples <- names(sort(covariances, 
                                decreasing=T)[1:min(numrefs, length(covariances))])

cat(paste('[',format(Sys.time()),']'," Calculating the mean RC ...",sep=''))
sample.mean.counts <- mean(counts[, case_sample])
sample.sumcounts <- apply(counts[, reference.samples], 2, sum)

# normalize reference samples to sample of interest
cat(paste('[',format(Sys.time()),']'," normalize reference samples to sample of interest\n",sep=''))
counts[, reference.samples] <- apply(counts[, reference.samples], 2, 
                                     function(x, sample.mean.counts) 
                                       round(x * sample.mean.counts / 
                                               mean(x)), sample.mean.counts)  

# select reference samples and weightings using non-negative least squares
cat(paste('[',format(Sys.time()),']'," select reference samples and weightings using non-negative least squares\n",sep=''))
b <- counts[, case_sample]
A <- as.matrix(counts[, reference.samples])
library(nnls)
all <- nnls(A, b)$x
est <- matrix(0, nrow=50, ncol=length(reference.samples))
set.seed(1)
for (i in 1:50){
  d <- sample(nrow(A), min(500, nrow(A)))
  est[i, ] <- nnls(A[d, ], b[d])$x
}
weights <- colMeans(est)
sample.weights <- weights / sum(weights)


w_matrix <- as.matrix(sample.weights)
ref_weighted_rc <- A%*%w_matrix

library(Hmisc)
# calculate weighted mean of read count
# this is used to calculate emission probabilities
cat(paste('[',format(Sys.time()),']'," calculate weighted mean of read count\n",sep=''))
canoes.reads$mean <- apply(canoes.reads[, reference.samples], 
                     1, wtd.mean, sample.weights)


targets <- canoes.reads$target
write.table(paste('[',format(Sys.time()),'] [',case_sample,'] exclude probes with all zero counts',sep=''), file=log_file, col.names=FALSE, append=TRUE)
cat(paste('[',format(Sys.time()),']'," exclude probes with all zero counts\n",sep=''))
# exclude probes with all zero counts
nonzero.rows <- canoes.reads$mean > 0
nonzero.rows.df <- data.frame(target=canoes.reads$target, 
                              nonzero.rows=nonzero.rows)

canoes.reads <- canoes.reads[nonzero.rows, ]

output_df <- cbind(RC_ref_weighted_df[,c(2,3,4)],ratio,RC_ref_weighted_df[,c(7,6,5)])
write.table(output_df, output_file_name, row.names=FALSE, col.names=FALSE,sep = '\t')


