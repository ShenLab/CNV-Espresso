## CANOES SPARK
cnv_espresso_transmission_analysis_file<- './cnv_info_w_img_tranmission_analysis_pilot_prediction.csv'
cnv_espresso_df = read.csv(cnv_espresso_transmission_analysis_file)

familyID_list <- unique(cnv_espresso_df$FamID)
print(paste("Number of family/trio:", length(familyID_list)))

transmission_rate_per_family <- function(cnv_espresso_df, familyID, q){
  cnv_espresso_df_select <- cnv_espresso_df[which(cnv_espresso_df$FamID == familyID),]
  num_parent_cnv_transmitted <- nrow(cnv_espresso_df_select[which(cnv_espresso_df_select$SQ>=q & cnv_espresso_df_select$Offspring_SQ>=q),])
  num_parent_cnv <- nrow(cnv_espresso_df_select[which(cnv_espresso_df_select$SQ>=q),])
  transmission_rate <- num_parent_cnv_transmitted/num_parent_cnv
  return(transmission_rate)
}

transmission_data = c()
for(q in seq(0,1,by=0.05)){
  print(q)
  for (familyID in familyID_list){
    transmission_rate <- transmission_rate_per_family(cnv_espresso_df, familyID, q)
    transmission_data <-rbind(transmission_data,data.frame(q, familyID, transmission_rate))
  }
}

data_path="."
pdf(paste(data_path,'/cnv_espresso_transmission_rate.pdf',sep=''))
boxplot(transmission_data$transmission_rate~transmission_data$q,boxwex=0.4, cex.lab=1.2, xlab = "Minimum SQ", ylab = "Transmission rate", main="")
abline(h=0.5, col = 'red')
dev.off()
cat("Transmission plot output to ",paste(data_path,'/cnvespresso_transmission_rate.pdf',sep=''))
