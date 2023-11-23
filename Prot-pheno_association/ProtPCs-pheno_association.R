setwd("/group/diangelantonio/users/alessia_mapelli/pQTL/Phenotype-Prot_associatioon")
library(FactoMineR)
library(lme4)
library(dplyr)
library(ggplot2)


file_path <- "/exchange/healthds/pQTL/INTERVAL/Proteomics_QC_files/data_all_not_imputed_without_transform_ANMLSMP_INTERVAL_QC.Rds"
data_new <- readRDS(file_path)
head(colnames(data_new), 33)
tail(colnames(data_new), 62)
col_of_start_protein <- 34
col_of_end_protein <- length(colnames(data_new)) - 61
protnames <- colnames(data_new)[col_of_start_protein:col_of_end_protein]
protdata <- data_new[protnames]
SampleId <- data_new$SampleId
export_prot_data <- cbind(protdata, SampleId)
write.csv(export_prot_data, "Interval_prot_data.csv")

phenonames_1 <- colnames(data_new)[1:col_of_start_protein-1]
phenodata_1 <- data_new[phenonames_1]
phenonames_2 <- colnames(data_new)[(col_of_end_protein+1): length(colnames(data_new))]
phenodata_2 <- data_new[phenonames_2]
phenodata <- cbind(phenodata_1, phenodata_2)
agePulse2 <- phenodata$agePulse*phenodata$agePulse
phenodata$bmi<-as.numeric(phenodata$wt_bl)/(as.numeric(phenodata$ht_bl)^2)
phenodata <- cbind(phenodata, agePulse2)
phenodata$processDate_bl <- as.Date(phenodata$processDate_bl, "%d%b%Y")
phenodata$attendanceDate <- as.Date(phenodata$attendanceDate, "%d%b%Y")
phenodata$ProcessSample <-
  as.POSIXct(paste(phenodata$processDate_bl, phenodata$processTime_bl, sep=" "), format = "%Y-%m-%d %H:%M:%S", tz="Europe/London")
phenodata$BloodDraw <-
  as.POSIXct(paste(phenodata$attendanceDate, phenodata$appointmentTime, sep=" "), format = "%Y-%m-%d %H:%M", tz="Europe/London")
phenodata$difftime <-
  as.numeric(difftime(phenodata$ProcessSample, phenodata$BloodDraw, units = "auto"))
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
write.csv(Pc, "PCs_10.csv")

Pc <- read.csv("PCs_10.csv")
Pc$SampleId <- as.factor(Pc$SampleId)
dtjoin<-left_join(Pc,phenodata,by="SampleId")
dev.off()
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.factor(Batch)))+geom_point()
p
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.factor(PlateId)))+geom_point()
p
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=ethnicPulse))+geom_point()
p
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.numeric(difftime)))+geom_point()
p
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.numeric(sexPulse)))+geom_point()
p
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.numeric(agePulse)))+geom_point()
p
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.numeric(bmi)))+geom_point()
p


dtjoin$sexPulse<-as.factor(dtjoin$sexPulse)
dtjoin$PlateId<-as.factor(dtjoin$PlateId)
dtjoin$Batch<-as.factor(dtjoin$Batch)
dtjoin$agePulse<-scale(dtjoin$agePulse)
dtjoin$bmi<-scale(dtjoin$bmi)
dtjoin$agePulse2<-scale(dtjoin$agePulse2)
dtjoin$difftime<-scale(dtjoin$difftime)
dtjoin$ethnicPulse<-as.factor(dtjoin$ethnicPulse)


mod<-lm(data=dtjoin,formula=Dim.1~difftime)
summary(mod)
mod<-lm(data=dtjoin,formula=Dim.2~difftime)
summary(mod)

mod<-lm(data=dtjoin,formula=Dim.1~sexPulse+agePulse+PlateId+difftime)
summary(mod)
mod<-lm(data=dtjoin,formula=Dim.2~sexPulse+agePulse+PlateId+difftime)
summary(mod)

