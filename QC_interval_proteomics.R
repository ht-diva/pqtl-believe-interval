# Required packages
packages <- c('usethis', 'devtools', 'SomaDataIO', 'FactoMineR', 'ggplot2', 'haven', 'dplyr', 'tidyverse', 'VCA', 'lme4')

# Install packages only if they are not already present
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)


library(haven)
library(SomaDataIO)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(FactoMineR)
library(VCA)
library(lme4)

################################################################################
################################################################################
###########USING ONLY SOMALOGIC DATA############################################

####Buffer, QC and min-max values extractions####
#loading proteomics data
#batch3
#my_adat <- read_adat("raw-data/SS-2227214_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
my_adat <-
  read_adat(
    "/processing_data/shared_datasets/plasma_proteome/believe/raw_data/proteome/batch_3/SS-2227214_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"
  )
is.soma_adat(my_adat)
metaprot<-getAnalyteInfo(my_adat)
metaprot$SeqId<-gsub("-",".", metaprot$SeqId)
metaprot$SeqId<-paste("seq.",metaprot$SeqId, sep="")
dim(my_adat)
##check samples notes
dim(my_adat[!is.na(my_adat$SampleNotes),])
##check if FAIL samples have sample notes
table(my_adat[!is.na(my_adat$SampleNotes),]$RowCheck)
b1<-as.data.frame(my_adat)
b1_all<-b1
dim(b1_all[b1_all$SampleType=="Sample",])
View(b1[!is.na(b1$SampleNotes),])
b1<-b1[b1$RowCheck=="PASS",]
dim(b1[b1$SampleType=="Sample",])
buffer_1<-b1[b1$SampleType=="Buffer",]
QC_1<-b1[b1$SampleType=="QC",]
calib_1<-b1[b1$SampleType=="Calibrator"&!is.na(b1$SampleType),]
list_target<-colnames(b1)[34:ncol(b1)]
data<-my_adat[,34:ncol(b1)]
dim(data)
table(is.na(data))
#table(is.numeric(data))
# summary(unlist(data))
# summary(unlist(b1[b1$SampleType=="Sample",34:ncol(b1)]))
#boxplot(log(unlist(data)))
#boxplot((unlist(data)))
# hist(log(unlist(data)))
##check max and min per target
max_1<-unlist(lapply(data, max))
min_1<-unlist(lapply(data, min))
plot(log(min_1),log(max_1))
which.max(max_1)
##create dataset with only numeric values for PASS samples
data<-data[rownames(b1),]
max_1<-unlist(lapply(data, max))
min_1<-unlist(lapply(data, min))
plot(log(min_1),log(max_1))
which.max(max_1)

rm(my_adat)
rm(data)
gc()

#batch4 same checks
# my_adat <- read_adat("raw-data/SS-2230701_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
my_adat <-
  read_adat(
    "/processing_data/shared_datasets/plasma_proteome/believe/raw_data/proteome/batch_4/SS-2230701_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"
  )
is.soma_adat(my_adat)
metaprot<-getAnalyteInfo(my_adat)
dim(my_adat)
dim(my_adat[!is.na(my_adat$SampleNotes),])
table(my_adat[!is.na(my_adat$SampleNotes),]$RowCheck)
b2<-as.data.frame(my_adat)
b2_all<-b2
dim(b2_all[b2_all$SampleType=="Sample",])
View(b2[!is.na(b2$SampleNotes),])
View(b2[is.na(b2$SampleType),])
b2<-b2[b2$RowCheck=="PASS"&!is.na(b2$SampleType),]
dim(b2[b2$SampleType=="Sample",])
buffer_2<-b2[b2$SampleType=="Buffer"&!is.na(b2$SampleType),]
QC_2<-b2[b2$SampleType=="QC"&!is.na(b2$SampleType),]
calib_2<-b2[b2$SampleType=="Calibrator"&!is.na(b2$SampleType),]
data<-my_adat[,34:ncol(b2)]
##ATTENTION targets start at column 35 for batch 2 (due to an additionnal non informative column)
dim(data)
table(is.na(data))
#table(is.numeric(data))
# summary(unlist(data))
# summary(unlist(b2[b2$SampleType=="Sample",34:ncol(b1)]))
# boxplot(log(unlist(data)))
# boxplot((unlist(data)))
# hist(log(unlist(data)))
max_2<-unlist(lapply(data, max))
min_2<-unlist(lapply(data, min))
plot(log(min_2),log(max_2))
which.max(max_2)
data<-data[rownames(b2),]
max_2<-unlist(lapply(data, max))
min_2<-unlist(lapply(data, min))
plot(log(min_2),log(max_2))
which.max(max_2)

rm(my_adat)
rm(data)
gc()


#####Mapping to genes#####
##count how many targets for the same gene
list_genes<-metaprot$EntrezGeneSymbol
list_genes[list_genes==""]<-NA
table(is.na(list_genes))
temp<-sort(table(list_genes))
temp<-temp%>% 
  as.data.frame() %>% 
  arrange(desc(Freq))
table(temp$Freq)

list_multiple<-list_genes[grep("\\|", list_genes)]

#merge QC and buffer to have one dataset for the 2 batches
QC_1$Batch<-rep(1,nrow(QC_1))
buffer_1$Batch<-rep(1,nrow(buffer_1))
calib_1$Batch<-rep(1,nrow(calib_1))

QC_2$Batch<-rep(2,nrow(QC_2))
buffer_2$Batch<-rep(2,nrow(buffer_2))
calib_2$Batch<-rep(2,nrow(calib_2))

buffer_all_row<-rbind(buffer_1,buffer_2)
QC_all_row<-rbind(QC_1,QC_2)
calib_all_row<-rbind(calib_1,calib_2)


####Min and max per batch####
max<-as.data.frame(rbind(max_1,max_2))
min<-as.data.frame(rbind(min_1,min_2))
MM<-as.data.frame(t(rbind(max, min)))
MM2<-MM
MM2$max<-pmax(MM2$max_1,MM2$max_2)
MM2$batch_max<-max.col(MM2[grepl("^max_", names(MM2))])
min.col <- function(m, ...) max.col(-m, ...)
MM2$min<-pmin(MM2$min_1,MM2$min_2)
MM2$batch_min<-min.col(MM2[grepl("^min_", names(MM2))])
MM2$SeqID<-row.names(MM2)

#plot
p<-ggplot(MM2,aes(x=log(min), y=log(max), color=as.factor(batch_max), ,shape=as.factor(batch_min)))+geom_point() +
  xlab("minimum log value for each protein")+ ylab("maxmum log value for each protein")+
  guides(color = guide_legend(title="Batch of the max value")) +
  guides(shape = guide_legend(title="Batch of the min value")) +xlim(-2,15)+ylim(-2,15)
p

####Batch effect on QC samples####

p.values <- lapply(34:7629, function(x,QC_all_row){
  # colnames(QC_all_row)[x]
  c(colnames(QC_all_row)[x],confint(lm(var1~batch,
                                       data=data.frame(cbind(var1=QC_all_row[,x],batch=as.factor(QC_all_row[,7630])))))[2,],
    summary(lm(var1~batch,
               data=data.frame(cbind(var1=QC_all_row[,x],batch=as.factor(QC_all_row[,7630])))))$coefficients[1,1],
    summary(lm(var1~batch,
               data=data.frame(cbind(var1=QC_all_row[,x],batch=as.factor(QC_all_row[,7630])))))$coefficients[2,])
},
QC_all_row)
p.values <- cbind(matrix(unlist(p.values), ncol = 8, byrow = TRUE)[,-7])
p.values<-as.data.frame(p.values)
colnames(p.values) <- c("var","conf - 2.5%","conf - 97.5%", "Intercept", "Est","Sd","pVal")
p.values <- p.values[p.values$var!="Intercept",]
p.values$pVal<-as.numeric(as.character(p.values$pVal))
p.values.adj<-p.values
pVal <- as.numeric(as.character(p.values$pVal))
wh <- which(p.adjust(pVal,"fdr")<0.05)
p.values.adj$pVal_adj<-p.adjust(pVal,"fdr")
wh <- p.values$var[wh]
length(wh)
p.values.adj$FC<-(as.numeric(p.values.adj$Intercept)+as.numeric(p.values.adj$Est))/as.numeric(p.values.adj$Intercept)
p.values.adj$Change<-abs(p.values.adj$FC-1)
p<-ggplot(p.values.adj,aes(y=var,x=-log(pVal_adj),col=(pVal<=0.05)))+geom_point()+ 
  geom_vline(xintercept = -log(0.05))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle("pValues adjusted for multiple comparisons - test of batch effect on QC samples")+
  ylab("Targets")
p
p<-ggplot(p.values.adj,aes(y=-log(pVal_adj),x=as.numeric(Est),col=(pVal<=0.05)))+geom_point()+ 
  geom_hline(yintercept = -log(0.05))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle("Volcano plot - test of batch effect on QC samples")+
  xlab("Effect estimates - log scale")
p
p<-ggplot(p.values.adj,aes(y=-log(pVal_adj),x=(as.numeric(FC)),col=(pVal<=0.05)))+geom_point()+ 
  geom_hline(yintercept = -log(0.05))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle("Volcano plot - test of batch effect on QC samples")+
  xlab("Fold change - log scale")
p

table(sign(as.numeric(p.values.adj$Est)))
table(sign(as.numeric(p.values.adj$Est[p.values.adj$pVal_adj<=0.05])))
p<-ggplot(QC_all_row,aes(x=seq.17145.1,y=seq.13378.80,col=as.factor(Batch)))+geom_point()
p
median(abs(as.numeric(p.values.adj$FC[p.values.adj$pVal_adj<=0.05])))
mean(abs(as.numeric(p.values.adj$FC[p.values.adj$pVal_adj<=0.05])))
median(p.values.adj$Change[p.values.adj$pVal_adj<=0.05])
mean(p.values.adj$Change[p.values.adj$pVal_adj<=0.05])
##compare directions of effect for target with highest significative change
table(p.values.adj$Change[p.values.adj$pVal_adj<=0.05]>0.05)
table(p.values.adj$Change[p.values.adj$pVal_adj<=0.05]>0.1)
table(p.values.adj$Change[p.values.adj$pVal_adj<=0.05]>0.2)
table(sign(as.numeric(p.values$Est[p.values.adj$Change>0.05&p.values.adj$pVal_adj<=0.05])))
table(sign(as.numeric(p.values$Est[p.values.adj$Change>0.1&p.values.adj$pVal_adj<=0.05])))
table(sign(as.numeric(p.values$Est[p.values.adj$Change>0.2&p.values.adj$pVal_adj<=0.05])))

table(as.numeric(p.values.adj$FC)>=1.05|as.numeric(p.values.adj$FC)<=0.95)
table(as.numeric(p.values.adj$FC)>=1.1|as.numeric(p.values.adj$FC)<=0.9)
table(as.numeric(p.values.adj$FC)>=1.2|as.numeric(p.values.adj$FC)<=0.80)
table(as.numeric(p.values.adj$FC)>=1.5|as.numeric(p.values.adj$FC)<=0.5)
table((as.numeric(p.values.adj$FC[p.values.adj$pVal_adj<=0.05])>=1.05|as.numeric(p.values.adj$FC[p.values.adj$pVal_adj<=0.05])<=0.95))
table(as.numeric(p.values.adj$FC)>=1.1|as.numeric(p.values.adj$FC)<=0.9)
table(as.numeric(p.values.adj$FC)>=1.2|as.numeric(p.values.adj$FC)<=0.80)
table(as.numeric(p.values.adj$FC)>=1.5|as.numeric(p.values.adj$FC)<=0.5)

##PCA on QC samples to assess batch effect
acp_p1<-PCA(QC_all_row[,34:(ncol(QC_all_row)-1)], ncp=10,scale.unit=TRUE,axes=c(1,2))
acp_p1$eig
QC_PC<-QC_all_row[,rownames(acp_p1$var$coord)]
A<-as.matrix(acp_p1$var$coord)
QC_PC<-as.matrix(QC_PC)
PC<-QC_PC%*%A
Pc<-as.data.frame(PC)
Pc$Batch<-QC_all_row$Batch
Pc$PlateID<-QC_all_row$PlateId
p<-ggplot(Pc,aes(x=Dim.1,y=Dim.2,col=as.factor(Batch)))+geom_point()
p
# p<-ggplot(Pc,aes(x=Dim.3,y=Dim.2,col=Batch))+geom_point()
# p
p<-ggplot(Pc,aes(x=Dim.4,y=Dim.5,col=PlateID))+geom_point()
p
# p<-ggplot(Pc,aes(x=Batch,y=Dim.8,col=Dim.3))+geom_point()
# p
p<-ggplot(Pc,aes(x=as.factor(Batch),y=Dim.8,col=Dim.3))+geom_boxplot()
p
mod<-lm(data=Pc,formula=Dim.1~Batch)
summary(mod)
mod<-lm(data=Pc,formula=Dim.2~Batch)
summary(mod)

####Computations of LOD#####
##see next section for plots
LOD_1<-lapply(buffer_1[34:(ncol(buffer_1)-1)], function(X)  {median(X)+5*mad(X)} )
LOD_2<-lapply(buffer_2[34:(ncol(buffer_2)-1)], function(X)  {median(X)+5*mad(X)} )
plot(log(as.numeric(LOD_1)),log(as.numeric(LOD_2)))
LOD<-as.data.frame(cbind(LOD1=unlist(LOD_1),LOD2=unlist(LOD_2)))
LOD$seq<-rownames(LOD)
View(LOD)     
LOD2<-pivot_longer(LOD,cols=c(1,2),names_to = "Batch", values_to = "LOD")
##correlation between LOD in the 2 batches
cor(LOD$LOD1,LOD$LOD2)

##Percent of values below LOD
b1<-b1[b1$SampleType=="Sample",]
table(unlist(lapply(1:(ncol(b1)-33), function(X) {
  colnames(b1)[X+33]==names(b1)[X]})))
underLOD1<-lapply(1:(ncol(b1)-33), function(X){
  table(b1[,33+X]<=LOD_1[[X]])["TRUE"]
})
names(underLOD1)<-names(LOD_1)
LOD_values<-cbind(LOD,unlist(underLOD1))
colnames(LOD_values)[4]<-"Samples_under_LOD_batch_1"
LOD_values$Percent_under_LOD_batch_1<-LOD_values$Samples_under_LOD_batch_1/nrow(b1)*100

b2<-b2[b2$SampleType=="Sample",]
table(unlist(lapply(1:(ncol(b2)-33), function(X) {
  colnames(b2)[X+33]==names(LOD_2)[X]})))
underLOD2<-lapply(1:(ncol(b2)-33), function(X){
  table(b2[,33+X]<=LOD_2[[X]])["TRUE"]
})
names(underLOD2)<-names(LOD_2)
table(rownames(LOD_values)==rownames(LOD))
LOD_values<-cbind(LOD_values,unlist(underLOD2))
colnames(LOD_values)[6]<-"Samples_under_LOD_batch_2"
LOD_values$Percent_under_LOD_batch_2<-LOD_values$Samples_under_LOD_batch_2/nrow(b2)*100
# saveRDS(LOD_values,"LOD_values_INTERVAL.Rds")
####CV buffer and min.max investigation in buffers####
CV_buffer<-lapply(1:(ncol(buffer_all_row)-34), function(X){
  temp<-as.numeric(buffer_all_row[,33+X])
  CV<-sd(temp)/mean(temp)
  return(CV)
})
CV_buffer_1<-lapply(1:(ncol(buffer_1)-34), function(X){
  temp<-as.numeric(buffer_1[,33+X])
  CV<-sd(temp)/mean(temp)
  return(CV)
})
CV_buffer_2<-lapply(1:(ncol(buffer_2)-34), function(X){
  temp<-as.numeric(buffer_2[,33+X])
  CV<-sd(temp)/mean(temp)
  return(CV)
})
max_buffer<-lapply(1:(ncol(buffer_all_row)-34), function(X){
  temp<-as.numeric(buffer_all_row[,33+X])
  return(max(temp))
})
index_max_buffer<-lapply(1:(ncol(buffer_all_row)-34), function(X){
  temp<-as.numeric(buffer_all_row[,33+X])
  return(which.max(temp))
})
min_buffer<-lapply(1:(ncol(buffer_all_row)-34), function(X){
  temp<-as.numeric(buffer_all_row[,33+X])
  return(min(temp))
})
index_min_buffer<-lapply(1:(ncol(buffer_all_row)-34), function(X){
  temp<-as.numeric(buffer_all_row[,33+X])
  return(which.min(temp))
})
plot(unlist(index_max_buffer),unlist(max_buffer))
plot(unlist(index_min_buffer),unlist(min_buffer))
names(CV_buffer)<-colnames(buffer_all_row)[34:(ncol(buffer_all_row)-1)]
table(names(CV_buffer)==rownames(LOD_values))
LOD_values$CV_from_buffer<-unlist(CV_buffer)
LOD_values$max_buffer<-unlist(max_buffer)
LOD_values$min_buffer<-unlist(min_buffer)
LOD_values_2<-LOD_values
LOD_values_2$index_max_buffer<-unlist(index_max_buffer)
LOD_values_2$check_of_sample_max_buffer<-buffer_all_row$RowCheck[unlist(index_max_buffer)]
LOD_values_2$index_min_buffer<-unlist(index_min_buffer)
LOD_values_2$check_of_sample_min_buffer<-buffer_all_row$RowCheck[unlist(index_min_buffer)]
p<-ggplot(LOD_values_2,aes(x=index_max_buffer,y=max_buffer,col=check_of_sample_max_buffer))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=index_max_buffer,y=max_buffer))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=index_max_buffer,y=log(max_buffer),col=check_of_sample_max_buffer))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=index_min_buffer,y=log(min_buffer),col=check_of_sample_min_buffer))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=index_min_buffer,y=(min_buffer)))+geom_point()
p
#each point corresponds to a target. Y is the maximun value assessed for this target in one of the blank. X is the identifier of the sample in which the maximum value was found
##for a specific target with an abnormal value  (same as in BELIEVE)
summary(buffer_all_row$seq.2171.12)
boxplot(buffer_all_row$seq.2171.12)
boxplot(log(buffer_all_row$seq.2780.35))
summary(unlist(buffer_all_row[,34:(ncol(buffer_all_row)-1)]))
boxplot(unlist(buffer_all_row[,34:(ncol(buffer_all_row)-1)]))
boxplot(log(unlist(buffer_all_row[,34:(ncol(buffer_all_row)-1)])))
##must also check with max seq.2171.12
#distribution LOD vs % under LOD
p<-ggplot(LOD_values_2,aes(x=LOD1,y=Percent_under_LOD_batch_1))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=LOD2,y=Percent_under_LOD_batch_2))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=log(LOD1),y=Percent_under_LOD_batch_1))+geom_point()
p
p<-ggplot(LOD_values_2,aes(y=seq,x=Percent_under_LOD_batch_1))+geom_point()+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  ylab("Targets")
p
p<-ggplot(LOD_values_2,aes(y=seq,x=Percent_under_LOD_batch_2))+geom_point()+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  ylab("Targets")
p

##Investagion of the patterns found in buffer and LOD according to the nature of the targets
###we expect the same signals in buffer and samples for the non-cleavable/hybridization control/non bioting and spuriomers aptamers
##selection of them
metaprot$SeqId<-gsub("-",".", metaprot$SeqId)
metaprot$SeqId<-paste("seq.",metaprot$SeqId, sep="")
list_human<-metaprot$SeqId[metaprot$Type=="Human"]
View(metaprot)
not_different<-metaprot$SeqId[metaprot$Type%in%c("Hybridization Control Elution","Non-Biotin","Spuriomer","Non-Cleavable")]
metaType<-metaprot[,c("Type","SeqId")]
LOD_values_2<-left_join(LOD_values_2,metaType,by=c("seq"="SeqId"))
p<-ggplot(LOD_values_2,aes(x=LOD1,y=Percent_under_LOD_batch_1,col=Type))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=LOD2,y=Percent_under_LOD_batch_2,col=Type))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=log(LOD1),y=Percent_under_LOD_batch_1,col=Type))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=log(LOD2),y=Percent_under_LOD_batch_2,col=Type))+geom_point()
p
p<-ggplot(LOD_values_2,aes(y=seq,x=Percent_under_LOD_batch_1,col=Type))+geom_point()+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  ylab("Targets")
p
p<-ggplot(LOD_values_2,aes(y=seq,x=Percent_under_LOD_batch_2,col=Type))+geom_point()+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  ylab("Targets")
p

p<-ggplot(LOD_values_2,aes(x=index_max_buffer,y=max_buffer,col=Type))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=index_max_buffer,y=max_buffer,col=Type))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=index_max_buffer,y=log(max_buffer),col=Type))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=index_min_buffer,y=log(min_buffer),col=Type))+geom_point()
p
p<-ggplot(LOD_values_2,aes(x=index_min_buffer,y=(min_buffer),col=Type))+geom_point()
p


table(LOD_values_2[LOD_values_2$Type=="Protein",]$Percent_under_LOD_batch_1>20)
table(LOD_values_2[LOD_values_2$Type=="Protein",]$Percent_under_LOD_batch_2>20)
table(LOD_values_2[LOD_values_2$Type=="Protein",]$Percent_under_LOD_batch_1>20|LOD_values_2[LOD_values_2$Type=="Protein",]$Percent_under_LOD_batch_2>20)


table(is.na(LOD_values_2[LOD_values_2$Type=="Protein",]$Percent_under_LOD_batch_1))
table(is.na(LOD_values_2[LOD_values_2$Type=="Protein",]$Percent_under_LOD_batch_2))
table(is.na(LOD_values_2[LOD_values_2$Type=="Protein",]$Percent_under_LOD_batch_1)|is.na(LOD_values_2[LOD_values_2$Type=="Protein",]$Percent_under_LOD_batch_2))

####CV QC####

#batch1
CV_QC_1<-lapply(1:(ncol(QC_1)-34), function(X){
  temp<-as.numeric(QC_1[,33+X])
  CV<-sd(temp)/mean(temp)
  return(CV)
})
table(unlist(CV_QC_1)<=0.2)
table(unlist(CV_QC_1)<=1)
#batch2
CV_QC_2<-lapply(1:(ncol(QC_2)-34), function(X){
  # a<-colnames(QC_2)[33+X]
  temp<-as.numeric(QC_2[,33+X])
  CV<-sd(temp)/mean(temp)
  # return(list(seq.id=a,CV=CV))
  return(CV)
})
table(unlist(CV_QC_2)<=0.2)
table(unlist(CV_QC_2)<=1)


CV_QC<-lapply(1:(ncol(QC_all_row)-34), function(X){
  temp<-as.numeric(QC_all_row[,33+X])
  CV<-sd(temp)/mean(temp)
  return(CV)
})
table(unlist(CV_QC)<=0.2)
table(unlist(CV_QC)<=1)

####CV calibrators###

#batch1
CV_calib_1<-lapply(1:(ncol(calib_1)-34), function(X){
  temp<-as.numeric(calib_1[,33+X])
  CV<-sd(temp)/mean(temp)
  return(CV)
})
table(unlist(CV_calib_1)<=0.2)
table(unlist(CV_calib_1)<=1)
#batch2
CV_calib_2<-lapply(1:(ncol(calib_2)-34), function(X){
  # a<-colnames(calib_2)[33+X]
  temp<-as.numeric(calib_2[,33+X])
  CV<-sd(temp)/mean(temp)
  # return(list(seq.id=a,CV=CV))
  return(CV)
})
table(unlist(CV_calib_2)<=0.2)
table(unlist(CV_calib_2)<=1)


CV_calib<-lapply(1:(ncol(calib_all_row)-34), function(X){
  temp<-as.numeric(calib_all_row[,33+X])
  CV<-sd(temp)/mean(temp)
  return(CV)
})
table(unlist(CV_calib)<=0.2)
table(unlist(CV_calib)<=1)

####CV calibrators intra plate####
##batch 1
plate_1<-calib_1$PlateId[!duplicated(calib_1$PlateId)]

CV_calib_1_intraplate<-lapply(1:length(plate_1), function(X){
  temp1<-calib_1[calib_1$PlateId==plate_1[X],]
  CV_X<-lapply(1:(ncol(temp1)-34), function(X){
    temp<-as.numeric(temp1[,33+X])
    CV<-sd(temp)/mean(temp)
    return(CV)
  })
  #return(list(plate_1[X],CV_X))
  return(CV_X)
})
CV_calib_1_intraplate<-as.data.frame(do.call("cbind", CV_calib_1_intraplate))
colnames(CV_calib_1_intraplate)<-plate_1
CV_calib_1_intraplate2<-as.data.frame(lapply(CV_calib_1_intraplate,unlist))
CV_calib1<-apply(CV_calib_1_intraplate2, 1, median, na.rm=T)
table(unlist(CV_calib1)<=0.2)
table(unlist(CV_calib1)<=1)

##batch2

plate_2<-calib_2$PlateId[!duplicated(calib_2$PlateId)]

CV_calib_2_intraplate<-lapply(1:length(plate_2), function(X){
  temp1<-calib_2[calib_2$PlateId==plate_2[X],]
  CV_X<-lapply(1:(ncol(temp1)-34), function(X){
    temp<-as.numeric(temp1[,33+X])
    CV<-sd(temp)/mean(temp)
    return(CV)
  })
  #return(list(plate_2[X],CV_X))
  return(CV_X)
})
CV_calib_2_intraplate<-as.data.frame(do.call("cbind", CV_calib_2_intraplate))
colnames(CV_calib_2_intraplate)<-plate_2
CV_calib_2_intraplate2<-as.data.frame(lapply(CV_calib_2_intraplate,unlist))
CV_calib2<-apply(CV_calib_2_intraplate2, 1, median, na.rm=T)
table(CV_calib2<=0.2)
table(CV_calib2<=1)

####Plot of CV and comparisons with Candia 2022####
CV_buf1<-as.numeric(unlist(CV_buffer_1))
CV_buf2<-as.numeric(unlist(CV_buffer_2))
CV_QC1<-unlist(CV_QC_1)
CV_QC2<-unlist(CV_QC_2)


##build one dataset with CV computed from buffer and QC  for all the targets
CV_all<-data.frame(cbind(targets=colnames(QC_1)[34:(ncol(QC_1)-1)], 
                         CV_buf1=as.numeric(CV_buf1), 
                         CV_buf2=as.numeric(CV_buf2), 
                         CV_QC1=as.numeric(CV_QC1),
                         CV_QC2=as.numeric(CV_QC2),
                         CV_calib1=as.numeric(CV_calib1),
                         CV_calib2=as.numeric(CV_calib2)))
CV_all$CV_buf1=as.numeric(CV_all$CV_buf1)
CV_all$CV_buf2=as.numeric(CV_all$CV_buf2)
CV_all$CV_QC1=as.numeric(CV_all$CV_QC1)
CV_all$CV_QC2=as.numeric(CV_all$CV_QC2)
CV_all$CV_calib1=as.numeric(CV_all$CV_calib1)
CV_all$CV_calib2=as.numeric(CV_all$CV_calib2)

##Compute buffer CV in Candia
# my_adat <- read_adat("raw-data/Candia/SS-217041.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
my_adat <- read_adat("/center/healthds/pQTL/Reference_datasets_for_QC_proteomics/Candia/SS-217041.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
View(my_adat)
table(my_adat$SampleType)
Candia_buffer<-as.data.frame(my_adat)
Candia_buffer<-Candia_buffer[Candia_buffer$SampleType=="Buffer",]
CV_buffer_candia<-lapply(1:(ncol(Candia_buffer)-33), function(X){
  temp<-as.numeric(Candia_buffer[,33+X])
  CV<-sd(temp)/mean(temp)
  return(list(colnames(Candia_buffer)[33+X],as.numeric(CV)))
})
mean_buffer_candia<-lapply(1:(ncol(Candia_buffer)-33), function(X){
  temp<-as.numeric(Candia_buffer[,33+X])
  m<-mean(temp)
  return(list(colnames(Candia_buffer)[33+X],as.numeric(m)))
})
mean_buffer_C<-do.call(rbind.data.frame, mean_buffer_candia)
colnames(mean_buffer_C)<-c("SeqId","mean_buffer_Candia")
CV_buffer_C<-do.call(rbind.data.frame, CV_buffer_candia)
colnames(CV_buffer_C)<-c("SeqId","CV_buffer_Candia")

##compute mean values in buffer in our data
mean_buffer_1<-lapply(1:(ncol(buffer_1)-34), function(X){
  temp<-as.numeric(buffer_1[,33+X])
  m<-mean(temp)
  return(list(colnames(buffer_1)[33+X],as.numeric(m)))
  
})
mean_buffer_b1<-do.call(rbind.data.frame, mean_buffer_1)
colnames(mean_buffer_b1)<-c("SeqId","mean_buffer_b1")
mean_buffer_2<-lapply(1:(ncol(buffer_2)-34), function(X){
  temp<-as.numeric(buffer_2[,33+X])
  m<-mean(temp)
  return(list(colnames(buffer_2)[33+X],as.numeric(m)))
  
})
mean_buffer_b2<-do.call(rbind.data.frame, mean_buffer_2)
colnames(mean_buffer_b2)<-c("SeqId","mean_buffer_b2")
mean_buffer<-left_join(mean_buffer_b1,mean_buffer_b2,by=c("SeqId"="SeqId"))
mean_buffer<-left_join(mean_buffer,mean_buffer_C,by=c("SeqId"="SeqId"))

##load CV computed from duplicates in Candia
# previousCV<-read.csv2("raw-data/SuppTab3Candia2022_without_header.csv",sep=",", header=T)
previousCV<-read.csv2("/center/healthds/pQTL/Reference_datasets_for_QC_proteomics/Candia/SuppTab3Candia2022_without_header.csv",sep=",", header=T)
previousCV$CV<-as.numeric(previousCV$hyb.msnCal.ps.cal.msnAll)
previousCV$SeqId<-gsub("-",".",previousCV$SeqId)
previousCV$SeqId<-paste("seq.",previousCV$SeqId,sep="")
##join dataset of CV from Candia and CV from our data
previousCV2<-left_join(previousCV,CV_buffer_C,by=c("SeqId"="SeqId"))
previousCV2$CV<-previousCV2$CV/100

dim(previousCV2)

##add the information about the dilution from the metadata from Somalogic
CV_all_2<-left_join(CV_all,previousCV2,by=c("targets"="SeqId"))
metaprot2<-metaprot[,colnames(metaprot)%in%c("SeqId","Dilution")]
# metaprot2$SeqId<-gsub("-",".", metaprot2$SeqId)
# metaprot2$SeqId<-paste("seq.",metaprot2$SeqId, sep="")


CV_all_2<-left_join(CV_all_2,metaprot2,by=c("targets"="SeqId"))
CV_all_2$Dilution.x<-as.numeric(CV_all_2$Dilution.x)
table(CV_all_2$Dilution.x==CV_all_2$Dilution.y)

##load information from Somalogic about expected QC
somalogicRef<-read.csv2("/center/healthds/pQTL/Reference_datasets_for_QC_proteomics/Somalogic/SomaScan_V4.1_7K_Annotated_Content_20210616.csv",sep=",", header=T)
somalogicRef$SeqId<-gsub("-",".",somalogicRef$SeqId)
somalogicRef$SeqId<-paste("seq.",somalogicRef$SeqId,sep="")
somalogicRef_IP<-somalogicRef[,c("SeqId","Intra.Plate.CV.Plasma")]
somalogicRef_IP$Intra.Plate.CV.Plasma<-as.numeric(somalogicRef_IP$Intra.Plate.CV.Plasma)
CV_all_2<-left_join(CV_all_2,somalogicRef_IP,by=c("targets"="SeqId"))
##various quantitative plots
p<-ggplot(CV_all_2,aes(x=log(CV_QC1),y=log(CV_QC2),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV computed on QC samples in the 2 batches")
p

p<-ggplot(CV_all_2,aes(x=log(CV_QC1),y=log(CV_calib1),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV computed on QC samples in batch 1 vs intraplate CV in batch 1")
p

p<-ggplot(CV_all_2,aes(x=log(CV_QC2),y=log(CV_calib2),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV computed on QC samples in batch 2 vs intraplate CV in batch 2")
p


p<-ggplot(CV_all_2,aes(x=log(CV_buf1),y=log(CV_calib1),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV computed on buffer samples in batch 1 vs intraplate CV in batch 1")
p

p<-ggplot(CV_all_2,aes(x=log(CV_buf2),y=log(CV_calib2),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV computed on buffer samples in batch 2 vs intraplate CV in batch 2")
p
p<-ggplot(CV_all_2,aes(x=log(CV_calib1),y=log(CV_calib2),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of  intraplate CV in the 2 batches")
p

p<-ggplot(CV_all_2,aes(x=log(CV_calib1),y=log(Intra.Plate.CV.Plasma),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of  intraplate CV in the batch 1 vs Somalogic references")
p

p<-ggplot(CV_all_2,aes(x=log(CV_calib2),y=log(Intra.Plate.CV.Plasma),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of  intraplate CV in the batch 2 vs Somalogic references")
p


p<-ggplot(CV_all_2,aes(x=log(CV_buf1),y=log(CV_buf2),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV computed on buffer samples in the 2 batches")

p

p1<-ggplot(CV_all_2,aes(y=log(CV_buf1),x=log(CV_buffer_Candia),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV in buffer - batch 1 vs Candia")
p1

p2<-ggplot(CV_all_2,aes(y=log(CV_buf2),x=log(CV_buffer_Candia),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV in buffer - batch 2 vs Candia")
p2
p3<-ggplot(CV_all_2,aes(x=log(CV),y=log(CV_QC1),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV in QC samples - batch 1 vs Candia")

p3
p4<-ggplot(CV_all_2,aes(x=log(CV),y=log(CV_QC2),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV in QC samples - batch 2 vs Candia")

p4


p<-ggplot(CV_all_2,aes(x=log(CV),y=log(CV_calib2),col=as.factor(Dilution.y)))+geom_point()+
  geom_hline(yintercept=log(0.2), linetype=2)+geom_vline(xintercept=log(0.2), linetype=2)+
  geom_abline(intercept=0, slope=1, linetype=1, color=2)+ggtitle("Comparisons of CV computedin Candia vs intraplate CV in batch 2")
p

table(CV_all_2$CV_QC1>=0.2)
table(CV_all_2$CV_QC2>=0.2)
table(CV_all_2$CV_buf1>=0.2)
table(CV_all_2$CV_buf2>=0.2)
table(CV_all_2$CV_QC1>=1)
table(CV_all_2$CV_QC2>=1)
table(CV_all_2$CV_buf1>=1)
table(CV_all_2$CV_buf2>=1)
#####cleaning workspace##########
#rm(QC_all_row)
rm(acp_p1)
rm(buffer_all_row)
rm(buffer_1)
rm(buffer_2)
rm(QC_1)
rm(QC_2)
rm(QC_PC)
rm(underLOD1)
rm(underLOD2)
rm(MM)
rm(min_nuffer)
rm(max_buffer)
rm(LOD_values)
rm(LOD_2)
#rm(LOD_values_2)
rm(LOD_1)
rm(LOD_2)
rm(LOD)
rm(index_max_buffer)
rm(index_min_buffer)

gc()

#####################################################################################
##Using phenotype covariates to compare biological and technical effects##
#####################################################################################

##load main covariates
# pheno<-read.csv("raw-data/Phenotype_interval/INTERVALdata_21DEC2022.csv")
pheno<-read.csv("/processing_data/shared_datasets/plasma_proteome/interval/phenotypes/INTERVALdata_21DEC2022.csv")
# link<-read.csv("raw-data/Phenotype_interval/INTERVAL_OmicsMap_20221221.csv")
link<-read.csv("/processing_data/shared_datasets/plasma_proteome/interval/phenotypes/INTERVAL_OmicsMap_20221221.csv")
table(duplicated(pheno$identifier))
table(duplicated(link$identifier))
all_pheno<-left_join(link,pheno,by="identifier")
table(!is.na(all_pheno$Soma7000_RAW))
pheno_g<-all_pheno[!is.na(all_pheno$Soma7000_RAW),]
table(duplicated(pheno_g$Soma7000_RAW))
pheno_g$bmi<-as.numeric(pheno_g$wt_bl)/(as.numeric(pheno_g$ht_bl)^2)
summary(pheno_g$ht_bl)
summary(pheno_g$wt_bl)
plot(pheno_g$ht_bl,pheno_g$wt_bl)
plot(pheno_g$ht_bl,pheno_g$agePulse)
plot(pheno_g$wt_bl,pheno_g$agePulse)
##pb with variable ht_bl, which is supposed to be self reported height in meter, but lots of very high and very low values and lots of identical values with lots of decimals
##check correspondence
#checking if Soma7000_raw from Phenotypes data corresponds to id in somalogic data
id_prot1<-b1$SampleId
id_prot1_WZ<-gsub("^0", "",id_prot1)
table(id_prot1%in%pheno_g$Soma7000_RAW)
table(id_prot1_WZ%in%pheno_g$Soma7000_RAW)
#all sampleID are in genid$4
table(duplicated(id_prot1_WZ))
#no duplicated gene ID

##prepare datasets
##batch1
dim(b1)
b1<-b1[b1$SampleType=="Sample",]
list_target<-colnames(b1)[34:ncol(b1)]
b1$id<-as.numeric(gsub("^0", "",as.character(b1$SampleId)))
dim(b1)
b1$batch<-rep(1,times=nrow(b1))
b1<-b1[,c("batch","id","PlateId","PlateRunDate","SlideId","SampleId","SampleType","SiteId","SubjectID","RowCheck",list_target)]
dim(b1)
##batch2
dim(b2)
# b2<-as.data.frame(my_adat)
b2<-b2[b2$SampleType=="Sample",]
# b2<-b2[b2$RowCheck=="PASS",]
# rm(my_adat)
b2$id<-as.numeric(gsub("^0", "",as.character(b2$SampleId)))
dim(b2)
b2$batch<-rep(2,times=nrow(b2))
b2<-b2[,c("batch","id","PlateId","PlateRunDate","SlideId","SampleId","SampleType","SiteId","SubjectID","RowCheck",list_target)]
dim(b2)
gc()

##create 1 big dataset
table(colnames(b1)==colnames(b2))
table(rownames(b1)%in%rownames(b2))
data<-rbind(b1,b2)
data_all<-left_join(data,pheno_g,by=c("id"="Soma7000_RAW"))

data_all<-as.data.frame(data_all)
data_all %>%count(batch,sexPulse)
data_all%>%count(batch,PlateId,sexPulse)
summary(data_all[data_all$batch==1,]$agePulse)
summary(data_all[data_all$batch==2,]$sexPulse)

####Sex Related Proteins#####
##using a list of sex related proteins provided by Somalogic
##sex related proteins
SRprot<-read.csv2("/center/healthds/pQTL/Reference_datasets_for_QC_proteomics/Somalogic/SexRelatedProteinsPlasma.csv",sep=",")
SRprot$seq.id<-gsub("-",".",SRprot$Sequence.ID)
SRprot$seq.id<-paste("seq.",SRprot$seq.id,sep="")
table(SRprot$seq.id%in%colnames(data_all))
data_SR<-data_all[,c("sexPulse", "batch","id","PlateId","PlateRunDate","SlideId","SampleId","SampleType","SiteId","SubjectID","RowCheck",SRprot$seq.id)]
p.values <- lapply(12:ncol(data_SR), function(x,data_SR){
  # colnames(data_SR)[x]
  c(colnames(data_SR)[x],confint(lm(var1~sexPulse,
                                    data=data.frame(cbind(var1=data_SR[,x],sexPulse=as.factor(data_SR$sexPulse)))))[2,],
    summary(lm(var1~sexPulse,
               data=data.frame(cbind(var1=data_SR[,x],sexPulse=as.factor(data_SR$sexPulse)))))$coefficients[2,])
},
data_SR)
p.values <- cbind(matrix(unlist(p.values), ncol = 7, byrow = TRUE)[,-6])
p.values<-as.data.frame(p.values)
colnames(p.values) <- c("var","conf - 2.5%","conf - 97.5%", "Est","Sd","pVal")
p.values <- p.values[p.values$var!="Intercept",]
p.values$pVal<-as.numeric(as.character(p.values$pVal))
p.values.adj<-p.values
pVal <- as.numeric(as.character(p.values$pVal))
wh <- which(p.adjust(pVal,"BH")<=0.05)
p.values.adj$pVal_adj<-p.adjust(pVal,"BH")
p.values.adj$pVal_bon<-p.values.adj$pVal*ncol(data)
wh <- p.values$var[wh]
length(wh)
p<-ggplot(p.values.adj,aes(y=var,x=-log(pVal_bon),col=(pVal<=0.05)))+geom_point()+ 
  geom_vline(xintercept = -log(0.05))+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle("pValues adjusted for multiple comparisons - test of sex effect on sex -related samples")+
  ylab("Targets")
p
p<-ggplot(data_SR,aes(x=seq.10743.13,y=seq.12077.32,col=as.factor(sexPulse)))+geom_point()
p
p<-ggplot(data_SR[data_SR$RowCheck=="PASS",],aes(x=seq.6580.29,y=seq.21232.39,col=as.factor(sexPulse)))+geom_point()+
  xlab("Pregnancy zone protein")+ylab("Benign Prostate specific Antigen")
p
#highest value for pregnacy zone protein for a male is for"380768467"
#highest value for pregnacy zone protein for a male is for"380760985"
p<-ggplot(data_SR,aes(x=seq.10743.13,y=as.factor(sexPulse)))+geom_boxplot()
p

##PCA for sex-related proteins
acp_p1<-PCA(data_SR[,12:ncol(data_SR)], ncp=10,scale.unit=TRUE,axes=c(1,2))
SR_PC<-data_SR[,rownames(acp_p1$var$coord)]
A<-as.matrix(acp_p1$var$coord)
SR_PC<-as.matrix(SR_PC)
PC<-SR_PC%*%A
Pc<-as.data.frame(PC)
Pc$sexPulse<-data_SR$sexPulse
p<-ggplot(Pc,aes(x=Dim.1,y=Dim.2,col=as.factor(sexPulse)))+geom_point()
p
data_SR$sexPulse<-as.factor(data_SR$sexPulse)
SRlonger<-data_SR %>% select(id,seq.6580.29,seq.21232.39,seq.8468.19,seq.7926.13,seq.13699.6,seq.4330.4,  
                             seq.3032.11,seq.16892.23,seq.2953.31,seq.9282.12,seq.5763.67,seq.7139.14,
                             seq.4914.10,seq.2474.54,seq.8428.102,seq.3396.54,seq.8484.24,seq.2575.5,
                             seq.3232.28,seq.5934.1,seq.15324.58,seq.20550.38,seq.5116.62,seq.10743.13,
                             seq.4874.3,seq.6973.111,seq.4234.8,seq.13744.37,seq.5735.54,seq.8989.40,
                             seq.3066.12,seq.12077.32,seq.3042.7,seq.13397.88,seq.8262.20,seq.6965.19,
                             seq.21112.6,seq.9002.36,seq.4929.55, sexPulse) %>%
  pivot_longer(., cols = c(seq.6580.29,seq.21232.39,seq.8468.19,seq.7926.13,seq.13699.6,seq.4330.4,  
                           seq.3032.11,seq.16892.23,seq.2953.31,seq.9282.12,seq.5763.67,seq.7139.14,
                           seq.4914.10,seq.2474.54,seq.8428.102,seq.3396.54,seq.8484.24,seq.2575.5,
                           seq.3232.28,seq.5934.1,seq.15324.58,seq.20550.38,seq.5116.62,seq.10743.13,
                           seq.4874.3,seq.6973.111,seq.4234.8,seq.13744.37,seq.5735.54,seq.8989.40,
                           seq.3066.12,seq.12077.32,seq.3042.7,seq.13397.88,seq.8262.20,seq.6965.19,
                           seq.21112.6,seq.9002.36,seq.4929.55), names_to = "Var", values_to = "Val")
for (i in (1:nrow(SRlonger))){
  SRlonger$label[i]<-SRprot$Protein.Name[SRprot$seq.id==SRlonger$Var[i]]
}
SRlonger$label<-SRprot$Protein.Name[SRprot$seq.id==SRlonger$Var]

p<- ggplot(SRlonger,aes(x = Var, y = log(Val), fill = sexPulse)) +
  geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
p<- ggplot(SRlonger,aes(x = label, y = log(Val), fill = sexPulse)) +
  geom_boxplot()+ theme(axis.text.x = element_text(angle = 45, hjust=1))
p

rm(SR_PC)
rm(SRlonger)
rm(SRprot)
gc()


####PCA to assess covariates effect####
##log transform, scale
# dim(data)
data[,11:ncol(data)]<-lapply(data[,11:ncol(data)],function(x) {scale(log(x))})
# dim(data)

#data<-data[data$batch=="1",]
#data<-data[data$batch=="2",]
acp_p1<-PCA(data[,11:ncol(data)], ncp=10,scale.unit=TRUE,axes=c(1,2))
dt_PC<-data[,rownames(acp_p1$var$coord)]
A<-as.matrix(acp_p1$var$coord)
dt_PC<-as.matrix(dt_PC)
PC<-dt_PC%*%A
table(rownames(PC)==rownames(data))
Pc<-as.data.frame(PC)
Pc<-cbind(Pc,data[,1:10])
dtjoin<-left_join(Pc,pheno_g,by=c("id"="Soma7000_RAW"))
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.factor(batch)))+geom_point()
p
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.factor(PlateId)))+geom_point()
p
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=ethnicPulse))+geom_point()
p
p<-ggplot(dtjoin,aes(x=Dim.1,y=Dim.2,col=as.numeric(bmi)))+geom_point()
p

dtjoin$sexPulse<-as.factor(dtjoin$sexPulse)
dtjoin$PlateId<-as.factor(dtjoin$PlateId)
dtjoin$batch<-as.factor(dtjoin$batch)
dtjoin$agePulse<-scale(dtjoin$agePulse)
dtjoin$bmi<-scale(dtjoin$bmi)
mod<-lm(data=dtjoin,formula=Dim.1~batch+sexPulse+agePulse)
summary(mod)
mod<-lm(data=dtjoin,formula=Dim.1~batch+sexPulse+agePulse+bmi+as.factor(smCurr_bl)+as.numeric(CRP_bl))
summary(mod)
mod<-lm(data=dtjoin,formula=Dim.1~batch+PlateId+difftime+sexPulse+agePulse+bmi+KidneyDisease)
summary(mod)
mod<-lm(data=dtjoin,formula=Dim.1~as.factor(PlateId)+difftime+sexPulse+agePulse+bmi+KidneyDisease)
summary(mod)

dtjoin$PlateId<-as.factor(dtjoin$PlateId)
mod<-lmer(Dim.1~batch+sexPulse+agePulse+(1|PlateId), data=dtjoin)
summary(mod)
mod<-lmer(Dim.1~batch+bmi+sexPulse+agePulse+as.factor(smCurr_bl)+as.numeric(CRP_bl)+(1|PlateId), data=dtjoin)
summary(mod)
# mod<-lmer(Dim.1~bmi+sexPulse+agePulse+KidneyDisease+(1|batch/PlateId), data=dtjoin)
# summary(mod)





