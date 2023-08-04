################################################################################
###Script to compare GWAS results for protein in condition A vs condition B ###
###############################################################################
library(dplyr)
library(tidyverse)
library(ggplot2)


####00. loading reference files#####

path<-"/center/healthds/pQTL/Solene_Believe_test/" ##to modify with your own path to file
#if the list of snp stratified by maf have not been created, create them from the 
##file of alternate  allele frequencies
##otherwise just load them
if (!file.exists(paste(path,"maf_005.txt",sep=""))){
  ##loading file with alternate allele frequencies
  filename<-"HDS_BELIEVE_final_HDS_maf_freq.afreq"
  af<-read.table(paste(path,filename,sep=""), header=F)
  af<-af[,c(2,5)]
  ##computing maf
  af$maf<-ifelse(af$V5<0.5,af$V5,1-af$V5)
  ##create 3 lists of snps corresponding to MAF>5%, 1%<=MAF<=5% and MAF<1%.
  maf_005<-af$V2[af$maf>0.05]
  maf_001<-af$V2[af$maf<=0.05&af$maf>=0.01]
  maf_min<-af$V2[af$maf<0.01]
  write
  rm(af)
  gc()
  ##saving results
  path<-"/center/healthds/pQTL/Solene_Believe_test/"
  write.table(maf_005,file=paste(path,"maf_005.txt",sep=""))
  write.table(maf_001,file=paste(path,"maf_001.txt",sep=""))
  write.table(maf_min,file=paste(path,"maf_min.txt",sep=""))
}else{
  ##loading maf results 
  maf_005<-(read.table(paste(path,"maf_005.txt",sep="")))$x
  maf_001<-(read.table(paste(path,"maf_001.txt",sep="")))$x
  maf_min<-(read.table(paste(path,"maf_min.txt",sep="")))$x
}



####0.loading GWAS results ####
#############
###regenie Cambridge###
#############

###############
####HT regenie pipe####
###############
##load top hits
#path to the folder of results
path<-"/group/diangelantonio/users/Solene/regenie_believe/test-pheno-replicate/" ##to modify with your own path to results
filenames <- list.files(paste(path,"results/tophits",sep=""), pattern="*.txt", full.names=TRUE) ##to modify with your own path to results
listOfFiles <- lapply(filenames, function(x) read.table(x, header = TRUE))
#1 file per residual tested
# names<-gsub('/group/diangelantonio/users/Solene/regenie_believe/test-pheno-replicate/results/tophits/', '', filenames)
names<-gsub(paste(path,'results/tophits/',sep=""), '', filenames)  ##to modify with your own path to results
names<-gsub('.regenie.filtered.annotated.txt', '', names) #to modify with your own name of results files
lapply(listOfFiles,dim)
names(listOfFiles)<-names
##function to retrieve the name of the model tested
identify_model<-function(listOfFiles){
  name<-names(listOfFiles)
  n<-length(listOfFiles)
  listUpdated<-list()
  for (i in (1:n)){
    temp<-listOfFiles[[i]]
    seq<-sub("\\_.*", "",name[i])
    model<-sub("^[^_]*_", "",name[i])
    temp$seq<-rep(seq,nrow(temp))
    temp$model<-rep(model,nrow(temp))
    listUpdated<-c(listUpdated,list(temp))
  }
  names(listUpdated)<-name
  return(listUpdated)
}
listOfFiles<-identify_model(listOfFiles)
df<-do.call("rbind",listOfFiles)
colnames(df)
table(df$model)##list of conditions for which we have hits reaching significant threshold

##choosing condition A
A<-"all_covar"
##choosing condition B
B<-"KidneyDisease"
###fastLMM###

##formatting results
##need to define the format

###load a MAF file for our genetic dataset
####1. inflation####

####2.GWAS signals####

####2.2. Venn diagram####
df_comp<-df[df$model%in%c(A,B),]
##targets for condition A
tar_A<-unique(df_comp$seq[df_comp$model==A])
##targets for condition B
tar_B<-unique(df_comp$seq[df_comp$model==B])
##all_targets
tar<-unique(c(tar_A,tar_B))
##TO ADD: check from the residuals files of targets to test if indeed all the targets from TAR were tested from the 2 proteins
##otherwise, remove them of the comparison

##Venn diagram by protein
snp_A<-lapply(tar_A, function(X)  {df_comp$ID[df_comp$model==A&df_comp$seq==X]})
names(snp_A)<-tar_A
snp_B<-lapply(tar_B, function(X)  {df_comp$ID[df_comp$model==B&df_comp$seq==X]})
names(snp_B)<-tar_B
common_hits_by_target<-function(X){
  if (X%in%tar_A&!X%in%tar_B){
    out<-c(length(snp_A[[X]]),0,0)
  }else{
    if (!X%in%tar_A&X%in%tar_B){
      out<-c(0,length(snp_B[[X]]),0)
    }else{
      a<-length(snp_A[[X]][snp_A[[X]]%in%snp_B[[X]]])
      out<-c(length(snp_A[[X]])-a,length(snp_B[[X]])-a,a)
    }
    
  }
  return(out)
}
list_Venn<-lapply(tar,common_hits_by_target)
df_Venn<-as.data.frame(do.call("rbind",list_Venn))
rownames(df_Venn)<-tar
df_Venn$seq<-rownames(df_Venn)
colnames(df_Venn)<-c("only_A","only_B","A_and_B","seq")
df_Venn$more_hits_in_model<-ifelse(df_Venn$only_A==df_Venn$only_B, "no_one",ifelse(df_Venn$only_A<df_Venn$only_B,"B","A"))
path_to_save<-"/center/healthds/pQTL/Solene_Believe_test/"
write.table(df_Venn,paste(path_to_save,"df_Venn_",A,"_",B,sep=""))


##Venn diagram across all protein stratified by cis and trans
##add here the attribution of cis and trans in df comp, once the file from Claudia is ready
##df_comp$cis_or_trans<-
df_Venn_stratified<-data.frame(categories=character(0),only_A=numeric(0),only_B=numeric(0), A_and_B=numeric(0))
##all targets
snp_A<-df_comp$ID[df_comp$model==A]
snp_B<-df_comp$ID[df_comp$model==B]
a<-length(snp_A[snp_A%in%snp_B])
df_Venn_stratified[1,]<-c("All",length(snp_A)-a,length(snp_B)-a,a)
# ##all targets cis
# a<-length(df_comp$ID[df_comp$model==A&df_comp$cis_or_trans=="cis"][df_comp$ID[df_comp$model==A&df_comp$cis_or_trans=="cis"]%in%df_comp$ID[df_comp$model==B&df_comp$cis_or_trans=="cis"]])
# df_Venn_stratified[2,]<-c("cis",length(df_comp$ID[df_comp$model==A&df_comp$cis_or_trans=="cis"])-a,length(df_comp$ID[df_comp$model==B&df_comp$cis_or_trans=="cis"])-a,a)
# ##all targets trans
# a<-length(df_comp$ID[df_comp$model==A&df_comp$cis_or_trans=="trans"][df_comp$ID[df_comp$model==A&df_comp$cis_or_trans=="trans"]%in%df_comp$ID[df_comp$model==B&df_comp$cis_or_trans=="trans"]])
# df_Venn_stratified[3,]<-c("trans",length(df_comp$ID[df_comp$model==A&df_comp$cis_or_trans=="trans"])-a,length(df_comp$ID[df_comp$model==B&df_comp$cis_or_trans=="trans"])-a,a)
##all targets maf_005
a<-length(snp_A[snp_A%in%maf_005&snp_A%in%snp_B])
df_Venn_stratified[4,]<-c("maf_005",length(snp_A[snp_A%in%maf_005])-a,length(snp_B[snp_B%in%maf_005])-a,a)
##all targets maf_001
a<-length(snp_A[snp_A%in%maf_001&snp_A%in%snp_B])
df_Venn_stratified[5,]<-c("maf_001",length(snp_A[snp_A%in%maf_001])-a,length(snp_B[snp_B%in%maf_001])-a,a)
##all_targets maf_min
a<-length(snp_A[snp_A%in%maf_min&snp_A%in%snp_B])
df_Venn_stratified[6,]<-c("maf_min",length(snp_A[snp_A%in%maf_min])-a,length(snp_B[snp_B%in%maf_min])-a,a)
df_Venn_stratified[,2:4]<-lapply(df_Venn_stratified[,2:4],as.numeric)
df_Venn_stratified$more_hits_in_model<-ifelse(df_Venn_stratified$only_A==df_Venn_stratified$only_B, "no_one",ifelse(df_Venn_stratified$only_A<df_Venn_stratified$only_B,"B","A"))
View(df_Venn_stratified)
#save results
path_to_save<-"/center/healthds/pQTL/Solene_Believe_test/"
write.table(df_Venn_stratified,paste(path_to_save,"df_Venn_stratified_",A,"_",B,sep=""))

####3. replication of previous studies####