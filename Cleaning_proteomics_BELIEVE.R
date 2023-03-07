#####Preparation dataset from raw data####

library(haven)
library(SomaDataIO)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(FactoMineR)
library(lme4)
library(labelled)
# compute distribution parameters for censored data
library(NADA)
# truncated normal distribution random draw (for fill-in method)
library(msm)
# intarval regression models
library(survival)
# wald test for survival regression terms
library(survey)
# to export survival model ps and betas as data.frame
library(broom)


################################
##Functions
protocol_survreg <- function(data_in) {
  # apply model and get p
  # defaults to log normal distribution
  fit1 <- survreg(
    Surv(val_start, val_end, type = "interval2") ~
      #Batch+difftime+sex+KidneyDisease+bmi,
      Batch,
    data = data_in,
    dist = "lognormal"
  )
  
  return(fit1)
}


replace_imputed <- function(original, imputed) {
  namestoChange <-
    colnames(imputed)[colnames(imputed) %in% colnames(original)]
  
  for (i in 1:length(namestoChange)) {
    original[namestoChange[i]] <- imputed[namestoChange[i]]
  }
  return(original)
  
}


################################
###load somalogic datasets
##batch 1
my_adat <-
  read_adat(
    "raw-data/SS-2218747_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat"
  )
is.soma_adat(my_adat)
metaprot1 <- as.data.frame(getAnalyteInfo(my_adat))
metaprot1 <- metaprot1[, c(
  "AptName",
  "SeqId",
  "SeqIdVersion",
  "SomaId",
  "TargetFullName",
  "Target",
  "UniProt",
  "EntrezGeneID",
  "EntrezGeneSymbol"  ,
  "Organism",
  "Units",
  "Type",
  "Dilution",
  "PlateScale_Reference",
  "CalReference",
  "ColCheck"
)]
b1 <- as.data.frame(my_adat)
list_target <- colnames(b1)[34:ncol(b1)]
list_covar <- colnames(b1)[1:33]
b1$Batch <- as.factor(rep(1, nrow(b1)))
list_covar <- c(list_covar, "Batch")
##keep only "PASS" samples
b1 <- b1[b1$RowCheck == "PASS",]
buffer_1 <- b1[b1$SampleType == "Buffer",]
QC_1 <- b1[b1$SampleType == "QC",]
calib_1 <-
  b1[b1$SampleType == "Calibrator" & !is.na(b1$SampleType),]
b1 <- b1[b1$SampleType == "Sample",]


##batch 2
my_adat <-
  read_adat(
    "raw-data/SS-2225082_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibration.anmlQC.qcCheck.anmlSMP.adat"
  )
is.soma_adat(my_adat)
my_adat <- my_adat[,-26]
metaprot2 <- as.data.frame(getAnalyteInfo(my_adat))
metaprot2 <- metaprot2[, c(
  "AptName",
  "SeqId",
  "SeqIdVersion",
  "SomaId",
  "TargetFullName",
  "Target",
  "UniProt",
  "EntrezGeneID",
  "EntrezGeneSymbol"  ,
  "Organism",
  "Units",
  "Type",
  "Dilution",
  "PlateScale_Reference",
  "CalReference",
  "ColCheck"
)]
b2 <- as.data.frame(my_adat)
b2$Batch <- as.factor(rep(2, nrow(b2)))
##keep only "PASS" samples
b2 <- b2[b2$RowCheck == "PASS" & !is.na(b2$SampleType),]
buffer_2 <- b2[b2$SampleType == "Buffer" & !is.na(b2$SampleType),]
QC_2 <- b2[b2$SampleType == "QC" & !is.na(b2$SampleType),]
calib_2 <-
  b2[b2$SampleType == "Calibrator" & !is.na(b2$SampleType),]
b2 <- b2[b2$SampleType == "Sample",]
rm(my_adat)
gc()


###merge -> main
table(colnames(buffer_1) == colnames(buffer_2))
table(colnames(QC_1) == colnames(QC_2))
table(colnames(calib_1) == colnames(calib_2))
table(colnames(b1) == colnames(b2))

buffer_all_row <- rbind(buffer_1, buffer_2)
QC_all_row <- rbind(QC_1, QC_2)
calib_all_row <- rbind(calib_1, calib_2)
main <- rbind(b1, b2)

table(duplicated(rownames(main)))


##build metadata dataframe
table(colnames(metaprot1) %in% colnames(metaprot2))
table(colnames(metaprot2) %in% colnames(metaprot1))
table(rownames(metaprot1) %in% rownames(metaprot2))

common_meta <-
  colnames(metaprot1)[(colnames(metaprot1) %in% colnames(metaprot2))]
all.equal(metaprot1[, common_meta], metaprot2[, common_meta])
names(metaprot1)[names(metaprot1) == 'ColCheck'] <- 'ColCheck1'
names(metaprot2)[names(metaprot2) == 'ColCheck'] <- 'ColCheck2'
metaprot <- merge(metaprot1, metaprot2, all.x = TRUE)
rm(metaprot1)
rm(metaprot2)
gc()

####rm nonhuman proteins from MD and main
table(metaprot$Organism)
##what about HIV-1 HIV-2? here keeping them
all_NH <- unique(metaprot$Organism)
all_NH <- all_NH[!all_NH %in% c("Human", "HIV-2", "HIV-2")]
all_NH_seq <- metaprot$AptName[metaprot$Organism %in% all_NH]
##updating list of target
list_target_updated <- list_target[!list_target %in% all_NH_seq]
length(list_target_updated)
##from MD
metaprot <- metaprot[!metaprot$AptName %in% all_NH_seq,]
dim(metaprot)
##from main
main <- main[,!colnames(main) %in% all_NH_seq]
dim(main)

##remove hybridization control, non cleavable, non biotin and spuriomer targets from MD and main
list_non_control_protein <-
  metaprot$AptName[metaprot$Type == "Protein"]
table(metaprot$Organism[metaprot$AptName %in% list_non_control_protein])
##updating list of target
list_target_updated <-
  list_target_updated[list_target_updated %in% list_non_control_protein]
##from MD
metaprot <-
  metaprot[metaprot$AptName %in% list_non_control_protein,]
dim(metaprot)
##from main
main <-
  main[, colnames(main) %in% c(list_covar, list_target_updated)]
dim(main)

###Compute LOD per batch
##formula from Candia 2022
LOD_1 <-
  lapply(list_target_updated, function(X) {
    median(buffer_1[, X]) + 5 * mad(buffer_1[, X])
  })
LOD_2 <-
  lapply(list_target_updated, function(X) {
    median(buffer_2[, X]) + 5 * mad(buffer_2[, X])
  })
LOD <-
  as.data.frame(cbind(LOD1 = unlist(LOD_1), LOD2 = unlist(LOD_2)))
LOD$AptName <- list_target_updated


##Percent of values below LOD
##batch 1
underLOD1 <- lapply(list_target_updated, function(X) {
  as.numeric(table(b1[, X] <= LOD$LOD1[LOD$AptName == X])["TRUE"]) / nrow(b1) *
    100
})
underLOD1 <-
  data.frame(Percent_under_LOD1 = unlist(underLOD1),
             AptName = list_target_updated)
##add FLAG to proteins with values below LOD in batch 1
underLOD1$FLAG_values_under_LOD_batch_1 <-
  !is.na(underLOD1$Percent_under_LOD1)
dim(underLOD1)

##batch 2
underLOD2 <- lapply(list_target_updated, function(X) {
  as.numeric(table(b2[, X] <= LOD$LOD2[LOD$AptName == X])["TRUE"]) / nrow(b2) *
    100
})
underLOD2 <-
  data.frame(Percent_under_LOD2 = unlist(underLOD2),
             AptName = list_target_updated)
##add FLAG to proteins with values below LOD in batch 1
underLOD2$FLAG_values_under_LOD_batch_2 <-
  !is.na(underLOD2$Percent_under_LOD2)
dim(underLOD2)

underLOD <- merge(underLOD1, underLOD2, all.x = T, all.y = T)
underLOD <- merge(underLOD, LOD, all.x = T)
dim(underLOD)
##remove targets with more than 50% of under LOD in at least one of the 2 batchs
##question: should we keep them and categorize them? or add an other FLAG?
underLOD <-
  underLOD[-which(underLOD$Percent_under_LOD1 >= 50 |
                    underLOD$Percent_under_LOD2 >= 50),]
list_target_updated <- underLOD$AptName
length(list_target_updated)
metaprot <- merge(underLOD, metaprot, all.x = T)
dim(metaprot)
main <-
  main[, colnames(main) %in% c(list_covar, list_target_updated)]
rm(underLOD)
rm(underLOD1)
rm(underLOD2)
gc()

###compute CV from QC

#batch1
CV_QC_1 <- lapply(list_target_updated, function(X) {
  temp <- as.numeric(QC_1[, X])
  CV <- sd(temp) / mean(temp) * 100
  return(CV)
})
CV_QC_1 <-
  data.frame(CV_QC1 = unlist(CV_QC_1), AptName = list_target_updated)
###add a variable flag QC in MD
CV_QC_1$FLAG_high_CV_1 <- (CV_QC_1$CV_QC1 >= 20)
#batch2
CV_QC_2 <- lapply(list_target_updated, function(X) {
  temp <- as.numeric(QC_2[, X])
  CV <- sd(temp) / mean(temp) * 100
  return(CV)
})
CV_QC_2 <-
  data.frame(CV_QC2 = unlist(CV_QC_2), AptName = list_target_updated)
###add a variable flag QC in MD
CV_QC_2$FLAG_high_CV_2 <- (CV_QC_2$CV_QC2 >= 20)
##merge
CV_QC <- merge(CV_QC_1, CV_QC_2, all.x = T, all.y = T)
metaprot <- merge(metaprot, CV_QC, all.x = T, all.y = T)
dim(metaprot)



##load variables from Believe

######Load BELIEVE covariates - prepare datasets############

##load main covariates
pheno <- read_dta("raw-data/Phenotypes_Cambridge/analysis.dta")
genid4 <- pheno$genid4
pheno_g <- pheno[!is.na(pheno$genid4),]
table(duplicated(pheno_g$genid4))
table(duplicated(pheno_g$idno))
pheno_g$genid2[pheno_g$genid2 == ""] <- NA
pheno_g$genid2[pheno_g$genid3 == ""] <- NA
pheno_g$genid2[pheno_g$genid == ""] <- NA
pheno_g$genid2[pheno_g$genid4 == ""] <- NA
table(is.na(pheno_g$genid))
table(is.na(pheno_g$genid2))
table(is.na(pheno_g$genid3))
table(is.na(pheno_g$genid4))
pheno_covar_BELIEVE <- colnames(pheno_g)

##load additionnal covariates
sup <-
  read.csv("raw-data/Phenotypes_Cambridge/BELIEVEdata_P5045_20221003_v2.csv")
sup$IDNO <- as.character(sup$IDNO)
sup$DateBloodDraw <- as.Date(sup$DateBloodDraw)
sup$DateLastMeal <- as.Date(sup$DateLastMeal)
sup$BloodDraw <-
  as.POSIXct(paste(sup$DateBloodDraw, sup$TimeBloodDraw), format = "%Y-%m-%d %H:%M:%S")
sup$LastMeal <-
  as.POSIXct(paste(sup$DateLastMeal, sup$TimeLastMeal), format = "%Y-%m-%d %H:%M:%S")
sup$KidneyDisease <- as.factor(sup$KidneyDisease)
##creation new variable difftime
sup$difftime <-
  as.numeric(difftime(sup$BloodDraw, sup$LastMeal, units = "mins"))
add_pheno_covar_BELIEVE <- colnames(sup)
pheno_covar_BELIEVE_updated <-
  c(pheno_covar_BELIEVE, add_pheno_covar_BELIEVE)
rm(pheno)
rm(buffer_1)
rm(buffer_2)
rm(calib_1)
rm(calib_2)
gc()
##merge
main$genid4 <-
  as.numeric(gsub("^0", "", as.character(main$SampleId)))

data_all <- left_join(main, pheno_g, by = c("genid4" = "genid4"))
data_all <- left_join(data_all, sup, by = c("idno" = "IDNO"))
data_all <- as.data.frame(data_all)
list_covar_updated <- c(list_covar, pheno_covar_BELIEVE_updated)
length(list_covar_updated)
list_covar_updated <-
  list_covar_updated[list_covar_updated %in% colnames(data_all)]
length(list_covar_updated)
data_all_not_imputed <- data_all
# saveRDS(data_all_not_imputed,"data_all_not_imputed_without_transform_Believe.Rds")



##imputation for variables below LOD
##Batch 1
data_for_imput <- data_all[data_all$Batch == 1, ]
data_for_imput <-
  data_for_imput[, c("SampleId", metaprot$AptName[metaprot$FLAG_values_under_LOD_batch_1])]

#from wide to long
imput_long <- data_for_imput %>%
  pivot_longer(cols = metaprot$AptName[metaprot$FLAG_values_under_LOD_batch_1],
               names_to = "target",
               values_to = "val")

# add LOD by joining dataset
imput_long <-
  left_join(imput_long, metaprot[, c("AptName", "LOD1")], by = c("target" =
                                                                   "AptName"))

# create concentration intervals for interval regression models (value start and
# end) and censorship indicators under lod yes/no
##attention il faut separer par batch
imput_long <- imput_long %>%
  mutate(
    val_start = case_when(val <= LOD1 ~ 0.1, ##pourquoi cette valeur
                          TRUE ~ val),
    val_end = case_when(val <= LOD1 ~ LOD1,
                        TRUE ~ val),
    under_LOD = val < LOD1
  )

# * fill-in ----
# compute distribution parameters for fill-in using cenros
##(regression on order statistics to model the distribution)
params <- imput_long %>%
  group_by(target) %>%
  do(stats_ros = NADA::cenros(log(.data$val), .data$under_LOD, forwardT = NULL)) %>%
  mutate(log_dist_mean = mean(stats_ros),
         log_dist_sd = sd(stats_ros)) %>%
  select(-stats_ros)



# add previously computed parameters to data
imput_long <- left_join(imput_long, params, by = "target")

# using random sample from truncated normal distribution named log_val_i for log
# transformed and imputed. data is drawn between log(val_start) and
# log(val_end)
##random sample generated from hidden Markov models (msm pacakge)
set.seed(22)
imput_long <- imput_long  %>%
  mutate(
    log_val_i = msm::rtnorm(
      n = nrow(imput_long),
      mean = log_dist_mean,
      sd = log_dist_sd,
      lower = log(val_start),
      upper = log(val_end)
    )
  )

imput_long <- imput_long %>%
  mutate(val_i = exp(log_val_i))
###pivot_wider
imput_final_1 <- imput_long %>%
  pivot_wider(
    id_cols = "SampleId",
    # id_cols = c("SampleId","sex","Batch","bmi","KidneyDisease","difftime"),
    names_from = c("target"),
    values_from = c("val_i")
  )


##Batch 2
data_for_imput <- data_all[data_all$Batch == 2, ]
data_for_imput <-
  data_for_imput[, c("SampleId", metaprot$AptName[metaprot$FLAG_values_under_LOD_batch_2])]

#from wide to long
imput_long <- data_for_imput %>%
  pivot_longer(cols = metaprot$AptName[metaprot$FLAG_values_under_LOD_batch_2],
               names_to = "target",
               values_to = "val")

# add LOD by joining dataset
imput_long <-
  left_join(imput_long, metaprot[, c("AptName", "LOD2")], by = c("target" =
                                                                   "AptName"))

# create concentration intervals for interval regression models (value start and
# end) and censorship indicators under lod yes/no

imput_long <- imput_long %>%
  mutate(
    val_start = case_when(val <= LOD2 ~ 0.1, ##pourquoi cette valeur
                          TRUE ~ val),
    val_end = case_when(val <= LOD2 ~ LOD2,
                        TRUE ~ val),
    under_LOD = val < LOD2
  )

# * fill-in ----
# compute distribution parameters for fill-in using cenros
##(regression on order statistics to model the distribution)
params <- imput_long %>%
  group_by(target) %>%
  do(stats_ros = NADA::cenros(log(.data$val), .data$under_LOD, forwardT = NULL)) %>%
  mutate(log_dist_mean = mean(stats_ros),
         log_dist_sd = sd(stats_ros)) %>%
  select(-stats_ros)

# add previously computed parameters to data
imput_long <- left_join(imput_long, params, by = "target")

# using random sample from truncated normal distribution named log_val_i for log
# transformed and imputed. data is drawn between log(val_start) and
# log(val_end)
##random sample generated from hidden Markov models (msm pacakge)
set.seed(22)
imput_long <- imput_long  %>%
  mutate(
    log_val_i = msm::rtnorm(
      n = nrow(imput_long),
      mean = log_dist_mean,
      sd = log_dist_sd,
      lower = log(val_start),
      upper = log(val_end)
    )
  )


imput_long <- imput_long %>%
  mutate(val_i = exp(log_val_i))



###pivot_wider
imput_final_2 <- imput_long %>%
  pivot_wider(
    id_cols = "SampleId",
    # id_cols = c("SampleId","sex","Batch","bmi","KidneyDisease","difftime"),
    names_from = c("target"),
    values_from = c("val_i")
  )


###merge keeping imputations
b1 <- data_all[data_all$Batch == 1, ]
table(b1$SampleId == imput_final_1$SampleId)
b1 <- replace_imputed(original = b1, imputed = imput_final_1)

b2 <- data_all[data_all$Batch == 2, ]
table(b2$SampleId == imput_final_2$SampleId)
b2 <- replace_imputed(original = b2, imputed = imput_final_2)

table(colnames(b1) == colnames(b2))
data_all_imputed <- rbind(b1, b2)


##reorder
data_all_imputed <- data_all_imputed %>%
  select(all_of(list_covar_updated), all_of(list_target_updated))
data_all_imputed_fillin_without_transform <- data_all_imputed
# saveRDS(data_all_imputed_fillin_without_transform,"data_all_imputed_fillin_without_transform_Believe.Rds")

#################################


###log transform all targets
###TO DO: SWITCH TO INVERSE NORMAL TRANSFORMATION?
#data_all_imputed[,(length(list_covar_updated)+1):ncol(data_all_imputed)]<-lapply(data_all_imputed[,(length(list_covar_updated)+1):ncol(data_all_imputed)],function(x) {scale(log(x))})
#data_all[,(length(list_covar_updated)+1):ncol(data_all)]<-lapply(data_all[,(length(list_covar_updated)+1):ncol(data_all)],function(x) {scale(log(x))})

####################################
####simple imputation of main covar
##prep median for correction covar
med_bmi <- median(data_all$bmi, na.rm = TRUE)
med_difftime <- median(data_all$difftime, na.rm = TRUE)

###prepare datasets for correction covar
data_all_imputed <- data_all_imputed  %>%
  mutate(
    sex = as_factor(sex),
    difftime = ifelse((difftime <= 0 |
                         is.na(difftime)), med_difftime, difftime),
    bmi = ifelse(is.na(bmi), med_bmi, bmi)
  )
data_all_imputed <- data_all_imputed %>%
  mutate_if(is.character, as.factor)

##save main and metaprot
#saveRDS(list(metadata=metaprot,imputed_cleaned_dataset=data_all_imputed),"QC_proteomics/cleaned_Believe.Rds")
saveRDS(
  list(metadata = metaprot, imputed_cleaned_dataset = data_all_imputed),
  "/center/healthds/pQTL/BELIEVE/cleaned_Believe.Rds"
)