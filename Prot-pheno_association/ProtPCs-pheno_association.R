setwd("/group/diangelantonio/users/alessia_mapelli/pQTL/Phenotype-Prot_associatioon")
library(FactoMineR)
library(lme4)
library(dplyr)
library(ggplot2)

col_of_start_protein <- 95

data_full <- readRDS("/center/healthds/pQTL/INTERVAL/cleaned_INTERVAL.Rds")$imputed_cleaned_dataset
protnames <- colnames(data_full)[col_of_start_protein:length(colnames(data_full))]
protdata <- data_full[protnames]
SampleId <- data_full$SampleId
export_prot_data <- cbind(protdata, SampleId)
write.csv(export_prot_data, "Interval_prot_data.csv")
phenonames <- colnames(data_full)[1:col_of_start_protein-1]
phenodata <- data_full[phenonames]
agePulse2 <- phenodata$agePulse*phenodata$agePulse
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
SampleId <- data_full$SampleId
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

mod<-lm(data=dtjoin,formula=Dim.1~Batch+sexPulse+agePulse+PlateId+bmi+difftime)
summary(mod)
mod<-lm(data=dtjoin,formula=Dim.2~Batch+sexPulse+agePulse+PlateId+bmi+difftime)
summary(mod)


mod<-lmer(Dim.1~Batch+sexPulse+agePulse+(1|PlateId), data=dtjoin)
summary(mod)
mod<-lmer(Dim.1~batch+bmi+sexPulse+agePulse+as.factor(smCurr_bl)+as.numeric(CRP_bl)+(1|PlateId), data=dtjoin)
summary(mod)
# mod<-lmer(Dim.1~bmi+sexPulse+agePulse+KidneyDisease+(1|batch/PlateId), data=dtjoin)
# summary(mod)
