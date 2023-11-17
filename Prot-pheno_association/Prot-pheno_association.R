library(survival)
library(dplyr)
library(ggrepel)
library(parallel)
install.packages(c("bigstatsr", "data.table"))
library(bigstatsr)
library(data.table)
library(ff)

################ UPLOAD THE EXPERIMENT CONFIGURATION FILE ###########

# set experiment config
configs <- RJSONIO::fromJSON("try_experiments_config.json")

conf = configs[[1]]

protdatapath <- conf$protdatapath
phenodatapath <- conf$phenodatapath
results_path <- conf$results_path
results_filename <- conf$results_filename
col_of_start_protein <- conf$col_of_start_protein
RDS <- conf$RDS

########### PREPARE THE DATASET

if(RDS == TRUE){
  dataset <- readRDS(protdatapath)$imputed_cleaned_dataset
  protnames <- colnames(dataset)[col_of_start_protein:length(colnames(dataset))]
  protdata <- dataset[protnames]
  phenonames <- colnames(dataset)[1:col_of_start_protein-1]
  phenodata <- dataset[phenonames]
  agePulse2 <- phenodata$agePulse*phenodata$agePulse
  phenodata <- cbind(phenodata, agePulse2)
  
} else{
  protdata <-read.csv(protdatapath)
  phenodata <-read.csv(phenodatapath)}


protdata_matrix <- as.matrix(protdata)
correlation_matrix <- cor(protdata_matrix)
save(correlation_matrix, file=paste(results_path,"/cor_matrix.Rdata", sep=""))
cor_matrix_rm <- correlation_matrix
cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0 # Modify correlation matrix
diag(cor_matrix_rm) <- 0

for(experiment_name in names(configs)){
  
  conf = configs[[experiment_name]]
  treatment_interaction = conf$treatment_interaction
  correlation_threshold = conf$correlation_threshold
  covariates = conf$covariates
  correction = conf$correction
  test_correction_n = conf$test_correction_n
  significance_thr = conf$significance_thr
  
  ############################
  covar_dataframe <- phenodata[covariates]
  final.ds = cbind(protdata, covar_dataframe)
  
  if(treatment_interaction == TRUE){
    summary_list <- mclapply(protnames, function(i) {
      system(paste("echo 'now processing: ",i,"'"))
      formula <- as.formula(sprintf("%s ~ %s + %s + %s*%s", i, covariates[1], covariates[2], covariates[1], covariates[2]))
      model <- lm(formula, data=final.ds, na.action=na.omit)
      data.frame(beta = summary(model)$coefficients[4,1],
                 pval = summary(model)$coefficients[4,4],
                 t_value = summary(model)$coefficients[4,3],
                 r_squared = as.numeric(summary(model)$r.squared))
    }
    )
    # construct one dataset out of lapply list of results
    result_ds = bind_rows(summary_list, .id = "column_label")[,-1]
    rownames(result_ds) <- prot_names
  }else{
      summary_list <- mclapply(protnames, function(i) {
      system(paste("echo 'now processing: ",i,"'"))
      formula <- as.formula(sprintf("%s ~ %s", i, covariates))
      model <- lm(formula, data=final.ds, na.action=na.omit)
      data.frame(beta = summary(model)$coefficients[2,1],
                 pval = summary(model)$coefficients[2,4],
                 t_value = summary(model)$coefficients[2,3],
                 r_squared = as.numeric(summary(model)$r.squared))
    }
    )
    # construct one dataset out of lapply list of results
    result_ds = bind_rows(summary_list, .id = "column_label")[,-1]
    rownames(result_ds) <- protnames
  }
  
  ### SIGNIFICANCE THRESHOLD (alpha) correcting for multiple testing #####
  
  # choose between effective number of tests, or the actual number of tests
  # the option test_correction_n = 'effective' estimates the effective number of tests based on Cheverud method [1]
  
  # set the number of tests in case we evaluate single patterns or patterns+patterns*ox_treatment
  if(treatment_interaction == TRUE){
    M = lenght(protnames)*3
  }else{
    M = length(protnames)
  }
  
  if(test_correction_n == 'effective'){
    corr_pca = prcomp(cor_matrix_rm)
    varlambdas = var((corr_pca$sdev)^2)
    M_eff = M*(1-(M-1)*varlambdas/M^2)
    N = M_eff
  } else {
    N = M
  }
  
  # chose the correction formula
  if(correction == 'bonferroni'){
    alpha = significance_thr/N
  } else {
    alpha = 1-(1-significance_thr)^(1/N)   ### Lynch more exact correction [2]
  }
  
  ###### check pvalue significance
  result_ds$signif_patt <- with(result_ds, ifelse(pval < alpha, TRUE, FALSE))
  significant_patterns = rownames(result_ds[result_ds$significant,])

  
  #### save results
  # create results directory if it doesn't exist yet
  if (!dir.exists(paste(results_path,"/results", sep=""))){
    dir.create(paste(results_path,"/results", sep=""))
  } 
  # save 
  result_file_name = paste(results_path,"/results/",results_filename,".RData", sep="")
  save(result_ds, significant_patterns, alpha, file = result_file_name)
  
}





### REFERENCES ####
# [1] Cheverud, J. M. A simple correction for multiple comparisons in interval mapping genome scans. Heredity. 87, 52–58 (2001).
# [2]  Lynch, M. and Walsh, B. (1998). Genetics and Analysis of Quantitative Traits. Sinauer Associates, Sunderland, MA.


