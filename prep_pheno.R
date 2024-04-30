############## This script reproduces parts of the pipeline to select covariates and inverse normalize the protein data, to create the phenotype file for Plink. 
############## TRIALS ON APRIL 16 2024
############## 
args <- commandArgs(trailingOnly=TRUE)
batch = as.character(args[1])

# using phenotypes from Xiyun and following steps from here to prepare phenotypes: /exchange/healthds/pQTL/results/INTERVAL/residuals/prepare_residuals.R
remove_ashk = FALSE
transform = TRUE
select_seqids = TRUE
residual = FALSE
split_by_batch = TRUE

if (select_seqids) {
	message("Provide a list of seqids to consider")
	# tested_seqid = c("seq.10708.3", "seq.4479.14")
	tested_seqids = c("seq.10361.25", "seq.10605.22", "seq.10708.3", "seq.11516.7", "seq.12618.50", "seq.14112.40", "seq.2780.35", "seq.2962.50", "seq.3815.14", "seq.3889.64", "seq.4479.14", "seq.5704.74", "seq.6493.9", "seq.7768.10", "seq.8068.43", "seq.8606.39", "seq.8877.22", "seq.8925.25", "seq.9398.30", "seq.9854.36")
}

if (split_by_batch) {
	message("Provide the batch to consider")
	#batch = 1
	print(batch)
}
###########################################
########## IMPORT PHENOTYPE
###########################################

### NOTE:
## This is a new file from Xiyun but we did not use this since the pQTLs were alsready completed
## data = readRDS("/exchange/healthds/pQTL/INTERVAL/residuals/INTERVAL_include_protein_LOD_proteomics/data_all_not_imputed_include_LOD_without_transform_ANMLSMP_INTERVAL_QC.Rds")

data = readRDS("/exchange/healthds/pQTL/INTERVAL/Proteomics_QC_files/data_all_not_imputed_without_transform_ANMLSMP_INTERVAL_QC.Rds")

###########################################
########## IMPORT GEN PC
###########################################

gen_PCs <- read.csv("/exchange/healthds/pQTL/INTERVAL/Genetic_QC_files/PCs_all_anc.csv")
gen_PCs <- gen_PCs[-1]

###########################################
########## SELECT AND TRANSFORM COVARIATES
###########################################

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

###########################################
########## SUBSET PHENO
###########################################
if (remove_ashk) {
	ashk = filtered_df_tot$FID[filtered_df_tot$imputed_anchestry=="Ashkenazi"]
	ashk_df = data.frame(FID=ashk, IID=ashk)
	write.table(ashk_df, file="/scratch/c.giambartolomei/TEST_INTERVAL_plink/april2024/152_individiauls_Ashkenazi.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}

covar = c("Batch", "agePulse", "sex", "difftime", "process_month", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
covardata = filtered_df_tot[, c("FID", covar)]


###########################################
########## SELECT PROTEIN DATA
###########################################

if (select_seqids) {
        protdata = filtered_df_tot[, c("FID", tested_seqids)]
}
if (!select_seqids) {
	all_seqids = names(filtered_df_tot)[grep("seq.", names(filtered_df_tot))]
        protdata = filtered_df_tot[, c("FID", all_seqids)]
}


###########################################
########## SELECT ALL OR SPLIT BY BATCH
###########################################

if (split_by_batch) {
	message("Selecting only batch ", batch)
	batch_fid = filtered_df_tot$FID[filtered_df_tot$Batch==batch]
	length(batch_fid)

	protdata = subset(protdata, protdata$FID %in% batch_fid)
	covardata = subset(covardata, covardata$FID %in% batch_fid)
}


###########################################
########## TRASNFROM PROTEIN DATA
###########################################

# This is from Michela/Solene/Alessia's script /exchange/healthds/pQTL/results/INTERVAL/residuals/protein_residuals_fitting.R
if (transform == TRUE){
## inverse transform the protein data
print("Data loaded.Transforming protein levels applying Rank-based Inverse Normal Transformation...")
# PROTEIN INV. RANK NORMAL TRANSFORM FUNCTION
inormal = function(x)
{
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

data_INV = protdata
data_INV = as.data.frame(lapply(data_INV[,-1], function(x) {scale(inormal(x))}))
data_INV = cbind.data.frame(FID = protdata[,"FID"], data_INV)
protdata = data_INV
}

protdata = cbind.data.frame(IID = protdata[,"FID"], protdata)

data = merge(protdata, covardata, by = "FID")

###########################################
########## If you want to continue creating the adjusted file used in the pQTL in Regenie 
########## Pipeline contiued from /exchange/healthds/pQTL/results/INTERVAL/residuals/protein_residuals_fitting.R
########## NOTE: with these specs, this file is identical to "/exchange/healthds/pQTL/results/INTERVAL/INTERVAL_NonImp_residuals_final.txt" that we used for pQTL analyses
#remove_ashk = FALSE
#transform = TRUE
#select_seqids = TRUE
#residual = TRUE
#split_by_batch = FALSE

if (residual) {

# REGRESSION FUNCTION
lm_reg_residuals = function(Y, dataset, model_definition)
  {
  reg=lm(paste(Y, model_definition), dataset, na.action = na.exclude)
  return(residuals(reg))
}

model_definition = "~ Batch + agePulse + sex + difftime + process_month + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
if (split_by_batch) {
	model_definition = "~ agePulse + sex + difftime + process_month + PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"
}

protnames = colnames(protdata)[-c(1,2)]
id_cols = c("IID", "FID")
data_all_info = data

  print("Starting with model fittings and residuals extraction.")

  counter = 0
  # Apply lm one protein at a time to obtain residuals with progress bar
  data_res <- lapply(protnames, function(protein) {
    # Increment the counter
    counter <<- counter + 1

    # Print a single line with the current progress
    cat("\rProcessing model", counter, "of", length(protnames), ". Target:", protein)

    lm_reg_residuals(Y = protein, dataset = data_all_info, model_definition = model_definition)
  })
  print("  ")
  print("Done with fitting linear models. Preparing final dataset and saving results.")

  ## put together the resulting dataframe
  data_res = as.data.frame(do.call("cbind", data_res))
  colnames(data_res)=protnames

  data_res = cbind(data_all_info[,id_cols], data_res)

  data = data_res
}


###########################################
########## WRITE
###########################################

#outfile_prefix = paste("pheno", "_remove_ashk", remove_ashk, "_transform", transform, "_residual", residual, sep="")
outfile_prefix = paste("pheno", "_transformINV_", transform, sep="")

if (split_by_batch) {
	outfile_prefix = paste(outfile_prefix, "_batch",  batch, sep = "")
}

if (residual) {
        outfile_prefix = paste(outfile_prefix, "_ResidualModel_", gsub(" |~", "", model_definition), sep="")
}

if (!residual) {
	outfile_prefix = paste(outfile_prefix, "_ResidualModel_FALSE", sep="")
}

if (remove_ashk) {
        outfile_prefix = paste(outfile_prefix, "remove_ashk", remove_ashk, sep="")
}

outfile = paste(outfile_prefix, ".txt", sep="")
message("Out file in: ", outfile)

write.table(data, file=outfile, sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

