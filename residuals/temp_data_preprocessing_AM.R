######################################################
if (!requireNamespace("imputeTS", quietly = TRUE)) {
  install.packages("imputeTS")
}
library(imputeTS)

# Impute NAs in the 'Value' column with the median

data = readRDS("/exchange/healthds/pQTL/INTERVAL/Proteomics_QC_files/data_all_not_imputed_without_transform_ANMLSMP_INTERVAL_QC.Rds")
gen_PCs <- read.csv("/exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/PCs_all_anc.csv")
gen_PCs <- gen_PCs[-1]

filtered_df <- data[data$Affymetrix_gwasQC_bl %in% gen_PCs$IID, ]

filtered_df$bmi<-as.numeric(filtered_df$wt_bl)/(as.numeric(filtered_df$ht_bl)^2)

filtered_df$processDate_bl <- as.Date(filtered_df$processDate_bl, "%d%b%Y")
filtered_df$attendanceDate <- as.Date(filtered_df$attendanceDate, "%d%b%Y")
filtered_df$ProcessSample <-
  as.POSIXct(paste(filtered_df$processDate_bl, filtered_df$processTime_bl, sep=" "), format = "%Y-%m-%d %H:%M:%S", tz="Europe/London")
filtered_df$BloodDraw <-
  as.POSIXct(paste(filtered_df$attendanceDate, filtered_df$appointmentTime, sep=" "), format = "%Y-%m-%d %H:%M", tz="Europe/London")
filtered_df$difftime <-
  as.numeric(difftime(filtered_df$ProcessSample, filtered_df$BloodDraw, units = "auto"))

filtered_df$difftime[is.na(filtered_df$difftime)] = median(filtered_df$difftime, na.rm = T)

filtered_df$sex = ifelse(filtered_df$sexPulse == 1, "M", "F")
filtered_df$process_month = format(filtered_df$processDate_bl, "%B")

filtered_df_tot = merge(filtered_df, gen_PCs, by.x = "Affymetrix_gwasQC_bl", by.y = "FID")
colnames(filtered_df_tot)[1] <- 'FID'

filtered_df_tot_controls <- filtered_df_tot [filtered_df_tot$SOMAPICK_CASE==0,]

write.table(filtered_df_tot, file = "INTERVAL_full_final.txt", sep='\t', row.names=FALSE, quote=FALSE, eol="\n")
write.table(filtered_df_tot_controls, file = "INTERVAL_contorl_final.txt", sep='\t', row.names=FALSE, quote=FALSE, eol="\n")


# gen_files = filtered_df_tot[, which(names(filtered_df_tot) %in% names(gen_PCs))]
# proteins_filtered = filtered_df[,34:7179]
# covariates = filtered_df[, -which(names(filtered_df) %in% names(proteins_filtered))]
# 
# write.table(covariates, file = "INTERVAL_covariates_final.txt", sep='\t', row.names=FALSE, quote=FALSE, eol="\n")
# write.table(proteins_filtered, file = "INTERVAL_proteins_NonImp_final.txt", sep='\t', row.names=FALSE, quote=FALSE, eol="\n")
# write.table(gen_PCs, file = "INTERVAL_geneticPCs_final.txt", sep='\t', row.names=FALSE, quote=FALSE, eol="\n")

# datapath <- "/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/INTERVAL_full_final.txt"
# data = read.table(datapath, header = TRUE, sep = '\t')
# colnames(data)[35]
# colnames(data)[7180]

