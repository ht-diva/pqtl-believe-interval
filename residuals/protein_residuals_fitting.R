############################################################
####### PROTEIN RESIDUALS FOR ASSOCIATION STUDIES ##########
############################################################


library(dplyr)
library(data.table)
library(tidyverse)
library(RJSONIO)



########## SET EXPERIMENT CONFIGS ##############
# NB: This is the only modification to make to this script.
# here one should define where the proper location of the config files.

conf_file_path = "residuals_config.json"

###############################################

configs = RJSONIO::fromJSON(conf_file_path)

# PROTEIN INV. RANK NORMAL TRANSFORM FUNCTION
inormal = function(x)
{
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# REGRESSION FUNCTION
lm_reg_residuals = function(Y, dataset, model_definition)
  {
  reg=lm(paste(Y, model_definition), dataset)
  return(reg$residuals)
}


number_of_configurations = length(configs)

for(experiment_name in names(configs)){

  ## parse configuration file and store all variables
  
  conf = configs[[experiment_name]]
  
  protdatapath = conf[['proteins_data_path']]
  covardatapath = conf[['covariates_data_path']]
  genpcpath = conf[['gen_PC_path']]
  datapath = conf[['unique_dataset_path']]
  results_path = conf[['results_path']]
  results_filename = conf[['results_filename']]
  n_start_prot = conf[['n_start_prot']]
  n_end_prot = conf[['n_end_prot']]
  n_datasets = conf[['n_datasets']]
  model_definition = conf[['model_definition']]
  add_suffix = conf[['add_suffix']]
  
  ## load data
  if (n_datasets > 1) {
    protdata = read.table(protdatapath, header = TRUE, sep = '\t')
    covariates = read.table(covardatapath, header = TRUE, sep = '\t')
    gen_PCs = read.table(genpcpath, header = TRUE, sep = '\t')
    
    ## collect names of proteins and samples in this experiment
    protnames = colnames(protdata)
    
  } else {
    data = read.table(datapath, header = TRUE, sep = '\t')
    protdata = data[,n_start_prot:n_end_prot]
    protnames = colnames(protdata)
    covariates = data[, -which(names(data) %in% names(protnames))]
  }
  
  ## inverse transform the protein data
  data_IVN = protdata
  data_IVN = as.data.frame(lapply(data_IVN, function(x) {scale(inormal(x))}))
  colnames(data_IVN)=protnames
  
  ## add all datasets together to give it as input to the linear regression models
  if (n_datasets > 1) {
    data_all_info = as.data.frame(cbind(data_IVN, gen_PCs, covariates))
  } else {
    data_all_info = as.data.frame(cbind(data_IVN, covariates))
  }
  
  # Apply lm one protein at a time to obtain residuals with progress bar
  data_res <- lapply(protnames, function(protein) {
    lm_reg_residuals(Y = protein, dataset = data_all_info, model_definition = model_definition)
  })

  ## put together the resulting dataframe
  data_res = as.data.frame(do.call("cbind", data_res))
  colnames(data_res)=protnames
  
  if (add_suffix == 1){
    colnames(data_res) = paste(colnames(data_res),"_res",sep="")
  }
  
  ## save resulting dataframe of residuals with the proper name according to the name
  ## of the configuration it is currently running
  if (number_of_configurations > 1) {
    results_filename = paste(results_filename, experiment_name, sep = "_")
  } else {
    results_filename = paste0(results_path, results_filename, '.txt')
  }
  
  write.table(data_res, file = results_filename, sep='\t', row.names=FALSE, quote=FALSE, eol="\n")
 
}


