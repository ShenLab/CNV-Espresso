#############################################################
# Author: Renjie Tan                                        #
# Date  : 3/2020                                            #
# This script is supplied as part of the CANOES pipeline.   #
# It has been modified to take in parameters from the shell # 
# script. It has been made more readable with additional    # 
# comments. It also sources the R script CANOES.R from the  #
# rscripts directory.                                       #
#############################################################
############################################################################
# STEP 1: Fetch the directory locations and file names from the shell script
############################################################################
library(doParallel)
library(foreach)
args <- commandArgs(trailingOnly=TRUE)
scripts_dir <- args[1]
data_dir <- args[2]
batch_name <- args[3]
output_dir <- args[4]
cores <- as.numeric(args[5])

#output_dir <- paste(output_dir, batch_name,'/',sep='')
output_file_name <- paste(output_dir,'canoes_calls_',batch_name,'.csv',sep='')#input_file_name <- paste(data_dir, '/cons_canoes_reads', sep='')
#sample_name_file <- paste(data_dir, '/sample_in_matrix.list', sep='')
#output_file_name <- paste(output_dir, '/canoes_calls_',batch_name,'.csv',sep='')

## For SPARK
input_file_name <- paste(data_dir, 'cons_canoes_reads_',batch_name,sep='') 
sample_name_file <- paste(data_dir,'merged_sample_list_', batch_name,'.list',sep='')
output_dir <- paste(output_dir, batch_name,'/',sep='')
output_file_name <- paste(output_dir,'canoes_calls_',batch_name,'.csv',sep='')

log_dir <- paste(output_dir,'/log/',sep='')
if(!dir.exists(output_dir)){
	dir.create(output_dir)
	dir.create(log_dir)
}

######################################################
# STEP 2: Read the read count file and GC content file
######################################################
setwd(data_dir)
gc <- read.table("gc.txt")$V2
canoes.reads <- read.table(input_file_name)

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

# RT moved the following scripts outside of CallCNVs function in CANOES.R.
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

## RT modified
# RT remove 'cov' variable to make sure the this variable is the right one for these samples.
# calculate covariance of read count across samples
#rm(covariance)
xcnv.list <- vector('list', length(sample.names))
if (file.exists(paste('./cov',basename(input_file_name),sep='.'))){
	load(file=paste('./cov',basename(input_file_name),sep='.'))
}
if (!exists('covariance')){
	cat(paste('[',format(Sys.time()),']'," calculate covariance of read count across samples\n",sep=''))
	system.time(covariance <- cor(canoes.reads[, sample.names], canoes.reads[, sample.names]))
	save(covariance,file=paste('./cov',basename(input_file_name),sep='.'))
}

######################################################
# STEP 4: Load R script and process each column/sample
######################################################
setwd(scripts_dir)
source("CANOES.R")
# cores <- detectCores(logical=F)
# cores <- floor(cores/10)
clus <- makeCluster(cores)
registerDoParallel(clus, cores=cores)
setwd(data_dir)
chunk.size <- length(sample.names)/cores
cat(paste("Start to run CANOES in ", as.character(cores), " cores ...\n",sep=''))
res_2.p <- foreach(i=1:cores, .combine='rbind') %dopar%
{
	log_file <- paste(log_dir,'/log_core_',i,'.log',sep='')
	write.table(paste('[',format(Sys.time()),'] cores:', i,' starts ....',sep=''), file=log_file, col.names=FALSE, append=FALSE)
	for (x in ((i-1)*chunk.size+1):(i*chunk.size)) {
        if (!file.exists(paste(output_dir,sample.names[x],sep=''))){
            cat(paste("Start to call CNVs on sample: ", sample.names[x], "\n",sep=''))
            xcnv.list[[x]] <- CallCNVs(sample.names[x], canoes.reads, covariance,log_file)
		    write.table(xcnv.list[[x]], file=paste(output_dir,sample.names[x],sep=''), col.names=FALSE, append=TRUE)
  		    gc()
        }    
	}
	xcnv.list
}
stopImplicitCluster()
stopCluster(clus)

#################################################################
# STEP 5: Append output as rows onto a single variable and write the output file
################################################################################
xcnvs <- do.call('rbind', xcnv.list)

setwd(data_dir)
#write.table(xcnvs, file=output_file_name, sep=",")

