install.packages("devtools")
devtools::install_github("SomaLogic/SomaDataIO")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("FactoMineR")
install.packages("VCA")
install.packages("haven")
install.packages("lme4")



install.packages("corrr")
library('corrr')
install.packages("ggcorrplot")
library(ggcorrplot)
library(haven)
library(SomaDataIO)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(FactoMineR)
library(VCA)
library(lme4)
library(resample)
library(mvtnorm)
library(RNOmni)

##################################
# Define the correct cohort
imputed<-readRDS("/center/healthds/pQTL/BELIEVE/cleaned_Believe.Rds")
data<-imputed$imputed_cleaned_dataset
###subset to ID with genetic QC data
fam_file<-read.table("/center/healthds/pQTL/BELIEVE/Genetic_QC_files/BELIEVE_genotype_final_bed_123456_withoutPCoutliers_mac.fam",sep="",header=F)
genid<-substr(fam_file$V1,19,27)
table(genid%in%data$genid)
rm(fam_file)
gc()
sum(!(data$genid%in%genid))
data<-data[data$genid%in%genid,]
saveRDS(data, file = "cleaned_Believe_restricted.rds")

######################################################
data<-read_csv("/group/diangelantonio/users/alessia_mapelli/PEER_prot_BELIEVE/cleaned_Believe_restricted.csv")

#######################################
# Compute the PCs

t_data <- apply(data[,68:ncol(data)],2, RankNorm)
hist(t_data[,10])
acp_p1<-PCA(t_data, ncp=20,scale.unit=TRUE,axes=c(1,2))
n_pca <-50
pov <- acp_p1$eig[1:n_pca,2]

#scree plot
barplot(pov, names.arg=1:n_pca, 
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")

plot(x = 1:n_pca, pov, names.arg=1:n_pca, 
     main = "Variances",
     xlab = "Principal Components",
     ylab = "Percentage of variances")
lines(x = 1:n_pca, pov, 
      type="b", pch=19, col = "steelblue")

cpv <- acp_p1$eig[1:n_pca,3]
plot(x = 1:n_pca, cpv, names.arg=1:n_pca, 
     main = "Variances",
     xlab = "Principal Components",
     ylab = "Percentage of variances")
lines(x = 1:n_pca, cpv, 
      type="b", pch=19, col = "steelblue")

###  
PC<-t_data[,rownames(acp_p1$var$coord)]
A<-as.matrix(acp_p1$var$coord)
PC<-as.matrix(PC)
PC<-PC%*%A
Pc<-as.data.frame(PC)
Pc$SampleId <- data$SampleId
write_csv(Pc, "PCs_20.csv")
Pc$Batch<-data$Batch
Pc$PlateID<-data$PlateId
dev.off()
p<-ggplot(Pc,aes(x=Dim.1,y=Dim.2))+geom_point()
p
p<-ggplot(Pc,aes(x=Dim.3,y=Dim.4))+geom_point()
p

mod<-lm(data=Pc,formula=Dim.1~Batch)
summary(mod)
mod<-lm(data=Pc,formula=Dim.2~Batch)
summary(mod)

summary(Pc)

