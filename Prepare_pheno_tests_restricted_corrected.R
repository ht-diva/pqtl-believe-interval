############################################################
###R script to prepare phenotypes files of targetsto test###
############################################################

library(dplyr)
library(data.table)
library(tidyverse)

################################################
###Load and prepare datasets with targets values
################################################

###load non imputed cleaned proteomic dataset
# non_imputed<-readRDS("/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/data_all_not_imputed_without_transform_Believe.Rds")
non_imputed<-readRDS("data_all_not_imputed_without_transform_Believe.Rds")
rownames(non_imputed)<-non_imputed$SampleId

###load imputed cleaned protoemics dataset
# imputed<-readRDS("/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/cleaned_Believe.Rds")
imputed<-readRDS("cleaned_Believe.Rds")
data<-imputed$imputed_cleaned_dataset
rownames(data)<-data$SampleId
table(colnames(data)%in%colnames(non_imputed))
non_imputed<-non_imputed[,colnames(data)]
table(colnames(data)==colnames(non_imputed))

##load non-ANML imputed cleaned dataset, restrict its number of samples/targets and reorder
# nANML<-readRDS("/group/diangelantonio/users/Solene/pQTL/Solene_Believe_test/cleaned_Believe_non_ANML.Rds")
nANML<-readRDS("cleaned_Believe_non_ANML.Rds")
nANML<-nANML$imputed_cleaned_dataset
rownames(nANML)<-nANML$SampleId
table(!colnames(nANML)%in%colnames(data))
table(!colnames(data)%in%colnames(nANML))
nANML<-nANML[,colnames(nANML)%in%colnames(data)]
dim(nANML)
dim(data)

nANML<-nANML[nANML$SampleId%in%data$SampleId,]
dim(nANML)
dim(data)
table(as.character(nANML$SampleId)==as.character(data$SampleId))
nANML<-nANML[as.character(data$SampleId),]
table(as.character(nANML$SampleId)==as.character(data$SampleId))

###subset all proteomics datasets to only ID with genetic QC data
##read genetic fam file
# fam_file<-read.table("/center/healthds/pQTL/Solene_Believe_test/BELIEVE_genotype_final_bed_123456_withoutPCoutliers_mac.fam",sep="",header=F)
fam_file<-read.table("BELIEVE_genotype_final_bed_123456_withoutPCoutliers_mac.fam",sep="",header=F)
genid<-substr(fam_file$V1,19,27)
table(genid%in%data$genid)
rm(fam_file)
gc()
non_imputed<-non_imputed[non_imputed$genid%in%genid,]
data<-data[data$genid%in%genid,]
nANML<-nANML[nANML$genid%in%genid,]
dim(data)
dim(non_imputed)
dim(nANML)
table(colnames(nANML)==colnames(data))
table(colnames(non_imputed)==colnames(data))


###Create 2 datasets with IVN and log transformation from imputed
inormal <- function(x)
{
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}
##variables 1 to 66 are covariates
data_IVN<-data
data_log<-data
rm(data)
rm(imputed)
gc()
non_imputed_IVN<-non_imputed
nANML_IVN<-nANML
data_IVN[,67:ncol(data_IVN)]<-lapply(data_IVN[,67:ncol(data_IVN)],function(x) {scale(inormal(x))})
data_log[,67:ncol(data_log)]<-lapply(data_log[,67:ncol(data_log)],function(x) {scale(log(x))})
non_imputed_IVN[,67:ncol(non_imputed_IVN)]<-lapply(non_imputed_IVN[,67:ncol(non_imputed_IVN)],function(x) {scale(inormal(x))})
nANML_IVN[,67:ncol(nANML_IVN)]<-lapply(nANML_IVN[,67:ncol(nANML_IVN)],function(x) {scale(inormal(x))})


#################################
##Prepare list of targets to test
#################################
##Prepare lists of cis and lowMaf from Claudia's list
# MergedDF = read.delim("/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/replicatedrstoseqid_sun_fenland_decode_add_chrom_pos.txt",header=T, sep="\t")
MergedDF = read.delim("replicatedrstoseqid_sun_fenland_decode_add_chrom_pos.txt",header=T, sep="\t")
table(MergedDF$SeqId.decode==MergedDF$SeqId.fenland)
table(MergedDF$SeqId.decode==MergedDF$SeqId.sun)
cis = MergedDF[MergedDF$cis == "cis",]
cis$SeqId.decode<-paste("seq.",cis$SeqId.decode,sep="")
table(cis$SeqId.decode%in%colnames(data_IVN))
cis$SeqId.decode[!cis$SeqId.decode%in%colnames(data_IVN)]
##7 targets are in the list from the literature but not our data -> to investigate
cis<-cis[cis$SeqId.decode%in%colnames(data_IVN),]

###Select 10 Low MAF targets
lowMaf<-cis[cis$MAF.sun<0.05 & cis$MAF.decode<0.05 & cis$MAF.fenland<0.05,c("ID", "SeqId.decode","Mapped.gene.fenland", "MAF.sun", "MAF.decode", "MAF.fenland")]
dim(lowMaf)
# 24 rsID-SeqIDs in cis have MAF lower than 0,05
lowMaf<-lowMaf[1:10,]
lowMaf$rsid<-gsub("_.*", "", lowMaf$ID)

###Select Top 10 cis 
cis_selected = cis[1:10,c("ID", "SeqId.decode","Mapped.gene.fenland", "MAF.sun", "MAF.decode", "MAF.fenland")]
cis_selected$rsid<-gsub("_.*", "", cis_selected$ID)# cis_selected[cis_selected$sun_literature,]
##check if no overlap between cis and low maf
table(cis_selected$rsid%in%lowMaf$rsid)
length(cis_selected[cis_selected$sun_literature,])
# write.table(lowMaf,file="/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/lowMaf.txt", sep='\t', row.names=FALSE,
#             #col.names=FALSE,
#             quote=FALSE, eol="\n")
# write.table(cis_selected,file="/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/cis_selected.txt", sep='\t', row.names=FALSE,
#             #col.names=FALSE,
#             quote=FALSE, eol="\n")
# write.table(as.data.frame(rbind(lowMaf,cis_selected)),file="/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/cis_maf.txt", sep='\t', row.names=FALSE,
#             #col.names=FALSE,
#             quote=FALSE, eol="\n")
write.table(as.data.frame(rbind(lowMaf,cis_selected)),file="cis_maf.txt", sep='\t', row.names=FALSE,
            #col.names=FALSE,
            quote=FALSE, eol="\n")

###Load list of trans from a list built manually
# trans_all<-read.csv2("/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/trans_replicated.csv",sep=",")
trans_all<-read.csv2("trans_replicated.csv",sep=",")
trans<-unique(trans_all$TargetId)
table(trans%in%colnames(data_IVN))

#############################################################
###Building a dataset with all targets to test and covariates
#############################################################
###Subsetting from the main datasets the targets we want to test 
##1.targets IVN with imputation/no imputation 
# list_targets<-read.csv2("/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/20230223_List_targets_to_test.csv",sep=",")
list_targets<-read.csv2("20230223_List_targets_to_test.csv",sep=",")
imputation_targets<-list_targets$SeqID[list_targets$Categories=="High % below LOD (10)"]
pheno<-data_IVN[,imputation_targets]
pheno_2<-non_imputed_IVN[,imputation_targets]
colnames(pheno_2)<-paste(colnames(pheno_2),"_non_imputed",sep="")
table(rownames(pheno)==rownames(pheno_2))
pheno<-as.data.frame(cbind(pheno,pheno_2))
##2. targets CV with IVN transformation
CV_targets<-list_targets$SeqID[list_targets$Categories=="High CV (10)"]
pheno_2<-data_IVN[,CV_targets]
table(rownames(pheno)==rownames(pheno_2))
pheno<-as.data.frame(cbind(pheno,pheno_2))
##3. for tests of models/transformation/anml
replicated_targets<-as.character(unique(c(lowMaf$SeqId.decode,cis_selected$SeqId.decode,trans)))
##attention 2 trans signals were already selected as lowMaf
##basic model on replicated targets (with values imputed, IVN transformed)
pheno_2<-data_IVN[,colnames(data_IVN)%in%replicated_targets]
rownames(pheno)
table(rownames(pheno)==rownames(pheno_2))
pheno<-as.data.frame(cbind(pheno,pheno_2))
##adding the non-ANML values for the replicated targets
pheno_2<-nANML_IVN[,colnames(nANML_IVN)%in%replicated_targets]
rownames(pheno)
table(rownames(pheno)==rownames(pheno_2))
colnames(pheno_2)<-paste(colnames(pheno_2),"_nANML",sep="")
pheno<-as.data.frame(cbind(pheno,pheno_2))
##adding the log transform values for the replicated targets
pheno_2<-data_log[,colnames(data_log)%in%replicated_targets]
table(rownames(pheno)==rownames(pheno_2))
colnames(pheno_2)<-paste(colnames(pheno_2),"_log",sep="")
pheno<-as.data.frame(cbind(pheno,pheno_2))

#### Adding all covar needed in the various models
basic_covar<-c("PlateId","Batch","ages","sex")
##question:shoudl we scale age?
optional_covar<-c("bmi","difftime","KidneyDisease","PlatePosition")
ID<-c("genid4")
##load genetic PCS
# gen_PCs<-read.table("/center/healthds/pQTL/Solene_Believe_test/BELIEVE_genotype_final_bed_PCs.txt", header=T)
gen_PCs<-read.table("BELIEVE_genotype_final_bed_PCs.txt", header=T)
gen_PCs$genid<-substr(gen_PCs$IID,19,27)
dim(gen_PCs)
table(gen_PCs$genid%in%data_IVN$genid)
rownames(gen_PCs)<-gen_PCs$genid
gen_PCs<-gen_PCs[as.character(data_IVN$genid),]
table(gen_PCs$genid==data_IVN$genid)
colnames(gen_PCs)
##add covar
pheno_2<-as.data.frame(cbind(data_IVN[,c(basic_covar,optional_covar,ID)],gen_PCs[,1:7]))
table(rownames(pheno)==rownames(pheno_2))
pheno<-as.data.frame(cbind(pheno,pheno_2))
##add proteomics PC
# Pc<-readRDS("/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/proteomics_PCs_common_subset.Rds")
Pc<-readRDS("proteomics_PCs_common_subset.Rds")
pheno_2<-as.data.frame(Pc)
#table(rownames(pheno)==rownames(pheno_2))
pheno<-as.data.frame(cbind(pheno,pheno_2))
##add PEER factors
# PEER<-read.table("/group/diangelantonio/users/alessia_mapelli/PEER_prot_BELIEVE/PEER 10/PEER_prot_BELIEVE_peers_10.txt",sep=",", header=TRUE)
PEER<-read.table("PEER_prot_BELIEVE_peers_10.txt",sep=",", header=TRUE)
PEER$genid4 <-
  as.numeric(gsub("^0", "", as.character(PEER$SampleId)))
rownames(PEER)<-PEER$genid4
pheno_2<-as.data.frame(PEER[,3:12])
table(pheno$genid4==rownames(pheno_2))
pheno<-as.data.frame(cbind(pheno,pheno_2))
#reorder
pheno_ordered<- pheno %>%
  dplyr::select(FID, IID, everything())
# write.table(pheno_ordered,file="/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/pheno_tests.txt", sep='\t', row.names=FALSE,
#             quote=FALSE, eol="\n")
View(pheno_ordered)

############################################################################
####Create a small explanatory with all targets to test and their categories
############################################################################
targets<-c(lowMaf$SeqId.decode,cis_selected$SeqId.decode,imputation_targets,CV_targets,trans)
length(c(lowMaf$SeqId.decode,cis_selected$SeqId.decode,imputation_targets,CV_targets,trans))
length(unique(c(lowMaf$SeqId.decode,cis_selected$SeqId.decode,imputation_targets,CV_targets,trans)))
table(trans%in%lowMaf$SeqId.decode)
##2 targets were both selected as trans and lowMAF
categories<-c(rep("lowMAF",length(lowMaf$SeqId.decode)),rep("cis",length(cis_selected$SeqId.decode)),
              rep("high%LOD", length(imputation_targets)),rep("highCV",length(CV_targets)),
              rep("trans",length(trans)))
df <- as.data.frame(cbind(targets=targets,categories=categories))
# write.csv2(df,"/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/Targets_to_test_categories_corrected.csv")
write.csv2(df,"Targets_to_test_categories_corrected.csv")




################################################################################
##preparation of a dataset directly with residuals of the various models to test
################################################################################
##recode binary in 0/1
pheno_cov<-pheno_ordered
pheno_cov$sex<-as.factor(ifelse(pheno_cov$sex=="Male",1,0))
pheno_cov$KidneyDisease<-as.factor(ifelse(pheno_cov$KidneyDisease=="2",0,1))
table(is.na(pheno_cov$KidneyDisease))
table(is.na(pheno_cov$bmi))
table(is.na(pheno_cov$difftime))
pheno_cov$Batch<-as.factor(ifelse(pheno_cov$Batch=="1",0,1))
##basic model (for all targets) ages sex batch plateId
reg<-function(X){
  reg1<-lm(paste(X,"~PlateId+ages+sex+Batch"),pheno_cov)
  return(reg1$residuals)
}
##remark: PlateId is nested inside Batch (removing BAtch wouldn't change anything)
# all_targets<-unique(targets) if we want to keep only the IVN imputed
all_targets<-colnames(pheno_cov)[3:116] #if we want to test non ANML/non imputation and log

pheno_res<-lapply(all_targets,reg)
pheno_res<-as.data.frame(do.call("cbind", pheno_res))
colnames(pheno_res)<-all_targets
colnames(pheno_res)<-paste(colnames(pheno_res),"_basic_covar",sep="")
pheno_3<-as.data.frame(cbind(FID=pheno_ordered$FID,IID=pheno_ordered$IID,pheno_res))
##add plate position only to replicated
reg<-function(X){
  reg1<-lm(paste(X,"~PlateId+ages+sex+Batch+PlatePosition"),pheno_cov)
  return(reg1$residuals)
}
pheno_res<-lapply(replicated_targets,reg)
pheno_res<-as.data.frame(do.call("cbind", pheno_res))
colnames(pheno_res)<-replicated_targets
colnames(pheno_res)<-paste(colnames(pheno_res),"_PlatePosition",sep="")
pheno_3<-as.data.frame(cbind(pheno_3,pheno_res))
##add difftime only to replicated
reg<-function(X){
  reg1<-lm(paste(X,"~PlateId+ages+sex+Batch+difftime"),pheno_cov)
  return(reg1$residuals)
}
pheno_res<-lapply(replicated_targets,reg)
pheno_res<-as.data.frame(do.call("cbind", pheno_res))
colnames(pheno_res)<-replicated_targets
colnames(pheno_res)<-paste(colnames(pheno_res),"_difftime",sep="")
pheno_3<-as.data.frame(cbind(pheno_3,pheno_res))
##add bmi only to replicated
reg<-function(X){
  reg1<-lm(paste(X,"~PlateId+ages+sex+Batch+bmi"),pheno_cov)
  return(reg1$residuals)
}
pheno_res<-lapply(replicated_targets,reg)
pheno_res<-as.data.frame(do.call("cbind", pheno_res))
colnames(pheno_res)<-replicated_targets
colnames(pheno_res)<-paste(colnames(pheno_res),"_bmi",sep="")
pheno_3<-as.data.frame(cbind(pheno_3,pheno_res))
##add KidneyDisease only to replicated
reg<-function(X){
  reg1<-lm(paste(X,"~PlateId+ages+sex+Batch+KidneyDisease"),pheno_cov)
  return(reg1$residuals)
}
pheno_res<-lapply(replicated_targets,reg)
pheno_res<-as.data.frame(do.call("cbind", pheno_res))
colnames(pheno_res)<-replicated_targets
colnames(pheno_res)<-paste(colnames(pheno_res),"_KidneyDisease",sep="")
pheno_3<-as.data.frame(cbind(pheno_3,pheno_res))
##add genetic PCs only to replicated
reg<-function(X){
  reg1<-lm(paste(X,"~PlateId+ages+sex+Batch+PC1+PC2+PC3+PC4+PC5"),pheno_cov)
  return(reg1$residuals)
}
pheno_res<-lapply(replicated_targets,reg)
pheno_res<-as.data.frame(do.call("cbind", pheno_res))
colnames(pheno_res)<-replicated_targets
colnames(pheno_res)<-paste(colnames(pheno_res),"_genPCs",sep="")
pheno_3<-as.data.frame(cbind(pheno_3,pheno_res))
##add protein PCs only to replicated
reg<-function(X){
  reg1<-lm(paste(X,"~PlateId+ages+sex+Batch+Dim.1+Dim.2+Dim.3+Dim.4+Dim.5+Dim.6+Dim.7"),pheno_cov)
  return(reg1$residuals)
}
pheno_res<-lapply(replicated_targets,reg)
pheno_res<-as.data.frame(do.call("cbind", pheno_res))
colnames(pheno_res)<-replicated_targets
colnames(pheno_res)<-paste(colnames(pheno_res),"_protPCs",sep="")
pheno_3<-as.data.frame(cbind(pheno_3,pheno_res))
###add PEER factors only to replicated
reg<-function(X){
  reg1<-lm(paste(X,"~PlateId+ages+sex+Batch+peer1+peer2+peer3+peer4+peer5+peer6+peer7+peer8+peer9+peer10"),pheno_cov)
  return(reg1$residuals)
}
pheno_res<-lapply(replicated_targets,reg)
pheno_res<-as.data.frame(do.call("cbind", pheno_res))
colnames(pheno_res)<-replicated_targets
colnames(pheno_res)<-paste(colnames(pheno_res),"_PEER",sep="")
pheno_3<-as.data.frame(cbind(pheno_3,pheno_res))
#model all covar (with prot PCs but not PEER)
reg<-function(X){
  reg1<-lm(paste(X,"~PlateId+ages+sex+Batch+PlatePosition+difftime+bmi+KidneyDisease+PC1+PC2+PC3+PC4+PC5+Dim.1+Dim.2+Dim.3+Dim.4+Dim.5+Dim.6+Dim.7"),pheno_cov)
  return(reg1$residuals)
}
pheno_res<-lapply(replicated_targets,reg)
pheno_res<-as.data.frame(do.call("cbind", pheno_res))
colnames(pheno_res)<-replicated_targets
colnames(pheno_res)<-paste(colnames(pheno_res),"_all_covar",sep="")
pheno_3<-as.data.frame(cbind(pheno_3,pheno_res))
#model all covar (with PEER but not protPCs)
reg<-function(X){
  reg1<-lm(paste(X,"~PlateId+ages+sex+Batch+PlatePosition+difftime+bmi+KidneyDisease+PC1+PC2+PC3+PC4+PC5+peer1+peer2+peer3+peer4+peer5+peer6+peer7+peer8+peer9+peer10"),pheno_cov)
  return(reg1$residuals)
}
pheno_res<-lapply(replicated_targets,reg)
pheno_res<-as.data.frame(do.call("cbind", pheno_res))
colnames(pheno_res)<-replicated_targets
colnames(pheno_res)<-paste(colnames(pheno_res),"_all_covar.PEER",sep="")
pheno_3<-as.data.frame(cbind(pheno_3,pheno_res))

dim(pheno_3)
colnames(pheno_3)

# write.table(pheno_3,file="/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/pheno_residuals_to_test.txt", sep='\t', row.names=FALSE,
#          quote=FALSE, eol="\n")
write.table(pheno_3,file="pheno_residuals_to_test_corrected.txt", sep='\t', row.names=FALSE,
            quote=FALSE, eol="\n")
# test<-read.table("/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/pheno_residuals_to_test.txt",header=T,sep='\t')
# test<-read.table("/home/solene.cadiou/QC_proteomics/pqtl-believe-interval/pheno_residuals_to_test_corrected.txt",header=T,sep='\t')

