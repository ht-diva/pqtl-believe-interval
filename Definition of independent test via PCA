setwd("/group/diangelantonio/users/alessia_mapelli/pQTL/INTERVAL/Phenotype-Prot_associatioon")
# traits is the phenotype input file for regenie 
traits <- read.csv("Interval_prot_data.csv")
# remove the index and SampleID columns
traits <- traits[,-c(1, dim(traits)[2])]

# base method
pc <- prcomp(traits, center=TRUE, scale.=TRUE )
summ <- summary(pc)
df <- summ$importance
idx <- which(df[3,] >= 0.95)[1]
cat("Method prcomp from base package determines", idx, "number of independent traits\n")
print(df[,idx])


# use the PCA method from the FactoMineR package
library(FactoMineR)
pca.imp <- PCA(traits)
df <- as.data.frame(pca.imp$eig)
idx <- which(df[,3] > 95)[1]
cat("Method PCR from FactoMinR package determins", idx, "number of independent traits\n")
print(df[idx,])

data<-as.data.frame(lapply(traits,function(x) {scale(x)}))
dim(data)

pca.imp <- PCA(data)
df <- as.data.frame(pca.imp$eig)
idx <- which(df[,3] > 95)[1]
cat("Method PCR from FactoMinR package with scaled data determins", idx, "number of independent traits\n")
print(df[idx,])

