setwd("/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Phenotype-Prot_associatioon")
library(FactoMineR)
library(lme4)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)


file_path <- "/exchange/healthds/pQTL/INTERVAL/Proteomics_QC_files/data_all_not_imputed_without_transform_ANMLSMP_INTERVAL_QC.Rds"
data_new <- readRDS(file_path)
dim(data_new)
protdata <- data_new %>%
            select(starts_with("seq."))
protnames <- colnames(protdata)
SampleId <- data_new$SampleId
export_prot_data <- cbind(protdata, SampleId)
write.csv(export_prot_data, "Interval_prot_data.csv")

phenodata <- data_new %>%
  select(-starts_with("seq."))

dim(phenodata)
phenodata$agePulse2 <- phenodata$agePulse*phenodata$agePulse
phenodata$bmi<-as.numeric(phenodata$wt_bl)/(as.numeric(phenodata$ht_bl)^2)
phenodata$processDate_bl <- as.Date(phenodata$processDate_bl, "%d%b%Y")
phenodata$attendanceDate <- as.Date(phenodata$attendanceDate, "%d%b%Y")
phenodata$ProcessSample <-
  as.POSIXct(paste(phenodata$processDate_bl, phenodata$processTime_bl, sep=" "), format = "%Y-%m-%d %H:%M:%S", tz="Europe/London")
phenodata$BloodDraw <-
  as.POSIXct(paste(phenodata$attendanceDate, phenodata$appointmentTime, sep=" "), format = "%Y-%m-%d %H:%M", tz="Europe/London")
phenodata$difftime <-
  as.numeric(difftime(phenodata$ProcessSample, phenodata$BloodDraw, units = "auto"))
phenodata$ProcessSample_month <- format(phenodata$ProcessSample, "%m")
phenodata$appointmentTime <- as.POSIXct(phenodata$appointmentTime, format = "%H:%M", tz="Europe/London")
phenodata <- phenodata %>%
  mutate(
  minutes_since_midnight = ifelse(
    is.na(appointmentTime), NA_integer_,
    as.integer(format(appointmentTime, "%H")) * 60 +
      as.integer(format(appointmentTime, "%M"))
  ),
  Time_of_day_blood_sample_collected = case_when(
    is.na(minutes_since_midnight)              ~ NA_character_,
    minutes_since_midnight <  8*60             ~ "early_morning",
    minutes_since_midnight < 12*60             ~ "morning",
    minutes_since_midnight < 16*60             ~ "afternoon",
    minutes_since_midnight < 19*60             ~ "late_afternoon",
    TRUE                                       ~ "night"
  ),
  Time_of_day_blood_sample_collected = factor(Time_of_day_blood_sample_collected,
                                              levels = c("early_morning","morning","afternoon",
                                                         "late_afternoon","night"),
                                              ordered = TRUE)
) %>%
  dplyr::select(-minutes_since_midnight)

write.csv(phenodata, "Interval_pheno_data.csv")
meta_data_full <- readRDS("/center/healthds/pQTL/INTERVAL/cleaned_INTERVAL.Rds")$metadata
write.csv(meta_data_full, "Interval_meta_data.csv")

####PCA to assess covariates effect####
##log transform, scale
data<-as.data.frame(lapply(protdata,function(x) {scale(log(x))}))
dim(data)

acp_p1<-PCA(data, ncp=10,scale.unit=TRUE,axes=c(1,2))
dt_PC<-data[,rownames(acp_p1$var$coord)]
A<-as.matrix(acp_p1$var$coord)
dt_PC<-as.matrix(dt_PC)
PC<-dt_PC%*%A
Pc<-as.data.frame(PC)
SampleId <- data_new$SampleId
Pc <- cbind(Pc, SampleId)
write.csv(Pc, "PCs_10_new.csv")

Pc <- read.csv("PCs_10_new.csv") %>%
  select( Dim.1, Dim.2, Dim.3, Dim.4, Dim.5, SampleId) %>%
  rename(Prot_PC1 =Dim.1,
         Prot_PC2=Dim.2, 
         Prot_PC3 = Dim.3,
         Prot_PC4=Dim.4,
         Prot_PC5=Dim.5)
pc_name = c("Prot_PC1","Prot_PC2","Prot_PC3","Prot_PC4","Prot_PC5")
phenodata <- read.csv("Interval_pheno_data.csv")
Pc$SampleId <- as.factor(Pc$SampleId)
phenodata$SampleId <- as.factor(phenodata$SampleId)
dtjoin<-left_join(Pc,phenodata,by="SampleId")

dtjoin$SampleId <- as.factor(dtjoin$SampleId)
dtjoin$sexPulse<-as.factor(dtjoin$sexPulse)
dtjoin$PlateId<-as.factor(dtjoin$PlateId)
dtjoin$Batch<-as.factor(dtjoin$Batch)
dtjoin$SOMAPICK_CASE<-as.factor(dtjoin$SOMAPICK_CASE)
dtjoin$agePulse<-scale(dtjoin$agePulse)
dtjoin$bmi<-scale(dtjoin$bmi)
dtjoin$agePulse2<-scale(dtjoin$agePulse2)
dtjoin$difftime<-scale(dtjoin$difftime)
dtjoin$ethnicPulse<-as.factor(dtjoin$ethnicPulse)
dtjoin$ProcessSample_month <- as.factor(dtjoin$ProcessSample_month)
dtjoin$smCurr_bl <- as.factor(dtjoin$smCurr_bl)
dtjoin$CRP_bl <- scale(dtjoin$CRP_bl)
dtjoin$Time_of_day_blood_sample_collected <- as.factor(dtjoin$Time_of_day_blood_sample_collected )


library(readr)
PCs_all_anc_new <- read_csv("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC/PCs_all_anc_new.csv")
colnames(PCs_all_anc_new)
View(dtjoin)

PCs_all_anc_new <- PCs_all_anc_new[,c("FID", "PC1", "PC2", "PC3", "PC4", "PC5",
                                      "PC6", "PC7", "PC8", "PC9", "PC10")]

colnames(PCs_all_anc_new) <- c("Affymetrix_gwasQC_bl","gen_PC1", "gen_PC2", "gen_PC3", "gen_PC4", "gen_PC5",
                               "gen_PC6", "gen_PC7", "gen_PC8", "gen_PC9", "gen_PC10" )

dtjoin <- left_join(dtjoin,PCs_all_anc_new,by="Affymetrix_gwasQC_bl")


dtjoin <- dtjoin %>%
  rename(Plate = PlateId,
         Time_btw_blood_sample_collected_and_sample_process = difftime, 
         Month_blood_sample_collected = ProcessSample_month,
         Soma_pick_case = SOMAPICK_CASE,
         Sex = sexPulse,
         Age = agePulse,
         Ethnicity = ethnicPulse,
         BMI = bmi,
         Smoking_current = smCurr_bl,
         CPR = CRP_bl,
         )

potential_conf_factors <- c(
  "Plate",
  "Batch",
  "Time_btw_blood_sample_collected_and_sample_process",
  "Month_blood_sample_collected",
  "Soma_pick_case",
  "Sex",
  "Age",
  "BMI",
  "Smoking_current",
  "CPR",
  "Time_of_day_blood_sample_collected",
  "gen_PC1", 
  "gen_PC2", 
  "gen_PC3", 
  "gen_PC4", 
  "gen_PC5",
  "KidneyDisease"
)

available_factors <- intersect(potential_conf_factors, colnames(dtjoin))
missing_factors <- setdiff(potential_conf_factors, colnames(dtjoin))

cat("Available confounding factors:", length(available_factors), "\n")
print(available_factors)

if (length(missing_factors) > 0) {
  cat("Missing confounding factors:", length(missing_factors), "\n")
  print(missing_factors)
}

cat("\n=== SYSTEMATIC CONFOUNDING FACTOR ANALYSIS ===\n")

# Function to calculate R-squared for PC ~ confounding factor
calculate_pc_rsquared <- function(pc_name, conf_factor, data) {
  tryCatch({
    # Handle different types of confounding factors
    if (is.character(data[[conf_factor]]) || is.factor(data[[conf_factor]])) {
      # For categorical variables
      formula_str <- paste(pc_name, "~", conf_factor)
    } else {
      # For continuous variables
      formula_str <- paste(pc_name, "~", conf_factor)
    }
    
    # Fit linear model
    model <- lm(as.formula(formula_str), data = data)
    model_summary <- summary(model)
    
    # Return adjusted R-squared (multiply by 100 for percentage)
    return(model_summary$adj.r.squared)
    
  }, error = function(e) {
    cat("Error calculating R² for", pc_name, "~", conf_factor, ":", e$message, "\n")
    return(NA)
  })
}

# Initialize results matrix
results_matrix <- matrix(NA, 
                         nrow = length(pc_name), 
                         ncol = length(available_factors),
                         dimnames = list(pc_name, available_factors))

cat("Calculating R² values for each PC-confounding factor combination...\n")

# Calculate R-squared for each combination
for (pc in pc_name) {
  for (factor in available_factors) {
    cat("Processing:", pc, "~", factor, "\n")
    rsq_value <- calculate_pc_rsquared(pc, factor, dtjoin)
    results_matrix[pc, factor] <- rsq_value
  }
}

# Convert to data frame for easier handling
results_df <- as.data.frame(results_matrix)
cat("Analysis completed!\n")


cat("\n=== CREATING SUMMARY HEATMAP ===\n")

# Prepare data for heatmap (transpose to match your example format)
heatmap_data <- t(results_matrix)*100

# Create a better color palette
# Use RdBu (Red-Blue) palette similar to your example
colors <- colorRampPalette(c( "white", "orange", "darkorange", "red"))(100)

# Create the main heatmap using pheatmap
pheatmap(
  heatmap_data,
  color = colors,
  breaks = seq(0, 
               max(abs(heatmap_data), na.rm = TRUE), 
               length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 8,
  main = "PC Confounding Factor Analysis\nR² Values (%)",
  fontsize = 10,
  angle_col = 45,
  cellwidth = 40,
  cellheight = 20,
  filename = "heatmap_plot_R_2_single_conf_factors.png",
  na_col = "grey90"
)

potential_multiple_conf_factors <- c(
  "Sex + Age + Time_btw_blood_sample_collected_and_sample_process + Month_blood_sample_collected",
  "Sex + Age + Time_btw_blood_sample_collected_and_sample_process + Month_blood_sample_collected + BMI",
  "Sex + Age + Time_btw_blood_sample_collected_and_sample_process + Month_blood_sample_collected + Plate",
  "Sex + Age + Time_btw_blood_sample_collected_and_sample_process + Month_blood_sample_collected + Batch",
  "Sex + Age + Time_btw_blood_sample_collected_and_sample_process + Month_blood_sample_collected + Plate + Batch",
  "Sex + Age + Time_btw_blood_sample_collected_and_sample_process + Month_blood_sample_collected + Batch + gen_PC1 + gen_PC2 + gen_PC3 + gen_PC4 + gen_PC5"
)

results_matrix_multiple <- matrix(NA, 
                         nrow = length(pc_name), 
                         ncol = length(potential_multiple_conf_factors),
                         dimnames = list(pc_name, potential_multiple_conf_factors))

cat("Calculating R² values for each PC-confounding factor combination...\n")

# Calculate R-squared for each combination
for (pc in pc_name) {
  for (factor in potential_multiple_conf_factors) {
    cat("Processing:", pc, "~", factor, "\n")
    rsq_value <- calculate_pc_rsquared(pc, factor, dtjoin)
    results_matrix_multiple[pc, factor] <- rsq_value
  }
}

results_df_mulitple <- as.data.frame(results_matrix_multiple)
cat("Analysis completed!\n")

View(t(results_df_mulitple))

results_df_full <- as.data.frame(rbind(t(results_df), t(results_df_mulitple)))
results_df_full <- results_df_full %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
View(results_df_full)
write.csv(results_df_full, "Non_genetic_assoc_ProtPCs.csv", row.names = TRUE)

install.packages("writexl") 
library(writexl)
writexl::write_xlsx(
  tibble::rownames_to_column(results_df_full, var = "Row"),
  path = "Non_genetic_assoc_ProtPCs.xlsx"
)

######### OLD CODE

# dev.off()
# p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.factor(Batch)))+geom_point()
# p
# p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.factor(PlateId)))+geom_point()
# p
# p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=ethnicPulse))+geom_point()
# p
# p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.numeric(difftime)))+geom_point()
# p
# p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.numeric(sexPulse)))+geom_point()
# p
# p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.numeric(agePulse)))+geom_point()
# p
# p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.numeric(bmi)))+geom_point()
# p
# p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.factor(SOMAPICK_CASE)))+geom_point()
# p
# p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.numeric(ProcessSample_month)))+geom_point()
# p

# #### IMPUTE NA with the MEDIAN
# colSums(is.na(dtjoin))
# dtjoin <- dtjoin[!is.na(dtjoin$ProcessSample_month), ] 
# dtjoin <- as.data.frame(lapply(dtjoin, function(x) ifelse(is.na(x), median(x, na.rm = TRUE), x)))
# colSums(is.na(dtjoin))
# write.csv(dtjoin, "PCs_pheno_no_nan.csv")
# 
# #### EXPLORATORY ANALYSIS
# hist(as.numeric(dtjoin$ProcessSample_month))
# hist(phenodata$difftime, xlab = "Difftime", main="Histogram of Difftime")
# 
# ### DEFINE THE MODELS
# # One for each single variable
# mod<-lm(data=dtjoin,formula=Dim.1~difftime)
# summary(mod)
# 
# # Multivaraite model
#mod<-lm(data=dtjoin,formula=Prot.PC4~Sex + Age + Time_btw_blood_sample_collected_and_sample_process + Month_blood_sample_collected)
#summary(mod)
# 
# 
# mod<-lm(data=dtjoin,formula=Dim.4~sexPulse+agePulse+ProcessSample_month+difftime+ Batch)
# summary(mod)
# 
# mod<-lm(data=dtjoin,formula=Dim.4~sexPulse+agePulse+ProcessSample_month+difftime+ Batch + gen.PC1+
#           gen.PC2+gen.PC3+gen.PC4+gen.PC5+gen.PC6+gen.PC7+gen.PC8+gen.PC9+gen.PC10)
# summary(mod)  

# ### DEFINE THE MODELS
# # One for each single variable
# mod<-lm(data=dtjoin,formula=Dim.1~difftime)
# summary(mod)
# 
# # Multivaraite model
#mod<-lm(data=dtjoin,formula=Prot.PC4~Sex + Age + Time_btw_blood_sample_collected_and_sample_process + Month_blood_sample_collected)
#summary(mod)
# 
# 
# mod<-lm(data=dtjoin,formula=Dim.4~sexPulse+agePulse+ProcessSample_month+difftime+ Batch)
# summary(mod)
# 

# Check Plate- Batch realtionship
dtjoin$Plate <- factor(dtjoin$Plate)
dtjoin$Batch <- factor(dtjoin$Batch)   # assumes it has values 1,2

# Contingency table
tab <- table(dtjoin$Plate, dtjoin$Batch)
tab
