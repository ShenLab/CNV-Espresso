## -----------------------------------------------------------------------------
#library(ggplot2)
library(plyr)
library(nnls)	

args <- commandArgs(trailingOnly=TRUE)
data_dir <- args[1]
test_sample_name <- args[2]
batch_name <- args[3]
output_dir <- args[4] 

numrefs <- 30
	
# sample list
sample_name_file <- paste(data_dir,'merged_sample_list_', batch_name,'.list',sep='')

# output file
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} 
output_file_name <- paste(output_dir,test_sample_name,'.txt',sep='')


## -----------------------------------------------------------------------------
# GC file
gc_file <- paste(data_dir, 'gc.txt.bz2',sep='')
if(file.exists(gc_file)){
    gc <- read.table(bzfile(gc_file))$V2
}else{
    stop("[Error] Could not find GC contend file.")
}

# input file
input_file_name <- paste(data_dir, 'cons_canoes_reads_',batch_name,'.bz2',sep='') 
input_file_RDate <- paste(input_file_name,'.RData',sep='') 
if(file.exists(input_file_RDate)){
  load(input_file_RDate)
}else{
  file.exists(input_file_name)
  canoes.reads <- read.table(bzfile(input_file_name))
  save(canoes.reads,file = input_file_RDate)
}


## -----------------------------------------------------------------------------
num_cols <- dim(canoes.reads)[2]
num_samples <- num_cols - 3
sample.names <- read.table(file=sample_name_file)
sample.names <- sub(".tar.bz2","",sample.names[,1])
sample.names <- sub(".cram","",sample.names)
sample.names <- sub(".bam","",sample.names)
names(canoes.reads) <- c("chromosome", "start", "end", sample.names)
target <- seq(1, nrow(canoes.reads))
canoes.reads <- cbind(target, canoes.reads[,1:3], gc, canoes.reads[,-c(1,2,3)])
#head(canoes.reads)


## -----------------------------------------------------------------------------
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
    canoes.reads$chromosome <- as.numeric(as.character(canoes.reads$chromosome))
}


## -----------------------------------------------------------------------------
# calculate covariance of read count across samples
file_name <- strsplit(basename(input_file_name),"\\.")[[1]][1]
cov_file  <- paste(data_dir,'cov.',file_name,sep='')
if (file.exists(cov_file)){
  load(file=cov_file)
}
if (!exists('covariance')){
  cat("Please calculate covariance of read count across samples at first.")
  cat(paste('[',format(Sys.time()),']'," Calculate covariance of read count across samples\n",sep=''))#
  system.time(covariance <- cor(canoes.reads[, sample.names], canoes.reads[, sample.names]))
  save(covariance, file=paste(data_dir, 'cov.', file_name, sep=''))
}


## -----------------------------------------------------------------------------
if (!test_sample_name %in% names(canoes.reads)){stop("No column for sample ", sample.name, " in counts matrix")}
if (length(setdiff(names(canoes.reads)[1:5], c("target", "chromosome", "start", "end", "gc"))) > 0){
  stop("First five columns of counts matrix must be target, chromosome, start, end, gc")
}
#nrow(canoes.reads[which(canoes.reads$chromosome=='NA'),])

if(class(canoes.reads$chromosome) != "numeric"){
  canoes.reads$chromosome <- as.numeric(as.character(canoes.reads$chromosome)) #RT modifided
}

canoes.reads <- arrange(canoes.reads, chromosome, start)

## -----------------------------------------------------------------------------
# find mean coverage of probes
cat(paste('[',format(Sys.time()),']'," Find mean coverage of probes\n",sep=''))
mean.counts <- mean(apply(canoes.reads[, sample.names], 2, mean))
#mean.counts

## -----------------------------------------------------------------------------
# normalize counts; round so we can use negative binomial
cat(paste('[',format(Sys.time()),']'," Normalize counts; round so we can use negative binomial\n",sep=''))
canoes.reads[, sample.names] <- apply(canoes.reads[, sample.names], 2, 
      function(x, mean.counts) 
               round(x * mean.counts / mean(x)), mean.counts)

# calculate covariance of read count across samples
reference.samples <- setdiff(sample.names, test_sample_name)
covariances <- covariance[test_sample_name, reference.samples]
reference.samples <- names(sort(covariances, 
        decreasing=T)[1:min(numrefs, length(covariances))])

## -----------------------------------------------------------------------------
cat(paste('[',format(Sys.time()),']'," Calculating the mean RC ...\n",sep=''))
sample.mean.counts <- mean(canoes.reads[, test_sample_name])
sample.sumcounts <- apply(canoes.reads[, reference.samples], 2, sum)

## -----------------------------------------------------------------------------
# normalize reference samples to sample of interest
cat(paste('[',format(Sys.time()),']'," normalize reference samples to sample of interest\n",sep=''))
canoes.reads[, reference.samples] <- apply(canoes.reads[, reference.samples], 2, 
      function(x, sample.mean.counts) 
              round(x * sample.mean.counts / 
              mean(x)), sample.mean.counts)  

## -----------------------------------------------------------------------------
# select reference samples and weightings using non-negative least squares
cat(paste('[',format(Sys.time()),']'," select reference samples and weightings using non-negative least squares\n",sep=''))
b <- canoes.reads[, test_sample_name]
A <- as.matrix(canoes.reads[, reference.samples])
all <- nnls(A, b)$x
est <- matrix(0, nrow=50, ncol=length(reference.samples))
set.seed(1)
for (i in 1:50){
  d <- sample(nrow(A), min(500, nrow(A)))
  est[i, ] <- nnls(A[d, ], b[d])$x
}
weights <- colMeans(est)
sample.weights <- weights / sum(weights)


## -----------------------------------------------------------------------------
library(Hmisc)
# calculate weighted mean of read count
cat(paste('[',format(Sys.time()),']'," Calculate weighted mean of read count\n",sep=''))
canoes.reads$weightedMean <- apply(canoes.reads[, reference.samples], 
                     1, wtd.mean, sample.weights) #RT equal to A%*%as.matrix(sample.weights)


## -----------------------------------------------------------------------------
cat(paste('[',format(Sys.time()),']'," Exclude probes with all zero counts\n",sep=''))
# exclude probes with all zero counts
nonzero.rows <- canoes.reads$weightedMean > 0
nonzero.rows.df <- data.frame(target=canoes.reads$target, 
                              nonzero.rows=nonzero.rows)
canoes.reads <- canoes.reads[nonzero.rows, ]

cat(paste('[',format(Sys.time()),']'," Calculate read count ratio.\n",sep=''))
canoes.reads$RCratio <- canoes.reads[,which(colnames(canoes.reads)==test_sample_name)] / canoes.reads$weightedMean
output_df <- cbind(canoes.reads[,c(2,3,4)],canoes.reads$RCratio,canoes.reads[,which(colnames(canoes.reads)==test_sample_name)],canoes.reads$weightedMean, canoes.reads$gc)
colnames(output_df) <- c("Chromosome","Start","End","RC_ratio", test_sample_name, "Weighted_Mean","GC")
write.table(output_df, output_file_name, row.names=FALSE, col.names=FALSE,sep = '\t')
cat(paste('[',format(Sys.time()),']'," Done.\n",sep=''))
