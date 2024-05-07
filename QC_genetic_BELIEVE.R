library(ggplot2)
library(haven)
library(dplyr)
library(tidyverse)
library(FactoMineR)
library(bigsnpr)
library(R.utils)
require(pgenlibr)
require(data.table)
library(bigsnpr)

devtools::source_url(
  "https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")
# path_to_plink<-"/home/solene.cadiou/plink2.sh"
path_to_plink<-"/plink2.sh"


###############################
####compute heterozigosity#####
###############################
het<-data.table::fread("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_subset_newHWE_het_12345.het",
                       sep = "\t",header = TRUE, verbose = TRUE)
View(het)
m <- mean(het$F) # Calculate the mean  
s <- sd(het$F) # Calculate the SD
valid <- (het$F <= m+3*s & het$F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
table(valid)
##excluse 284 individuals
hist(het$F,100)
abline(v=m+3*s)
abline(v=m-3*s)
to_remove_het<-as.data.frame(cbind(FID=het$`#FID`[!valid],IID=het$IID[!valid]))
write.table(to_remove_het,file="/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/GT_BELIEVE_subset_to_remove_het_HWE15.txt", sep='\t', row.names=FALSE,
            #col.names=FALSE,
            quote=FALSE, eol="\n")
##284

########################
#####SEX-CHECK##########
#########################
data <- read.table("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/bed_genotype_BELIEVE_subset_splitx_check_sex_fdistri_13_test_HWE15.sexcheck",sep="",header=T)
hist(data$F,100)
data$genid<-as.numeric(substr(data$IID,19,27))
data_2<-left_join(data,pheno, by = c("genid" = "genid"))

table((data_2$F<=0.35)&(data_2$sex==2))
table((data_2$F>=0.35)&(data_2$sex==1))
table(((data_2$F<=0.5)==(data_2$sex==2))&((data_2$F>=0.6)==(data_2$sex==1)))

obj.bed.sex<-snp_readBed("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/bed_genotype_BELIEVE_subset_splitx_imputesex_0506_13_test_HWE15.bed")
obj.bed.sex<-snp_attach("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/bed_genotype_BELIEVE_subset_splitx_imputesex_0506_13_test_HWE15.rds")
View(obj.bed.sex$fam)
sexcheck<-obj.bed.sex$fam
sexcheck$genid<-as.numeric(substr(sexcheck$sample.ID,19,27))
table(sexcheck$sex)
pheno <- read_dta("/processing_data/shared_datasets/plasma_proteome/believe/raw_data/phenotypes/analysis.dta")
pheno$genid[pheno$genid == ""] <- NA
sexcheck2<-left_join(sexcheck,pheno, by = c("genid" = "genid"))
table(sexcheck2$sex.x==sexcheck2$sex.y)
View(sexcheck2[!sexcheck2$sex.x==sexcheck2$sex.y,])
table(sexcheck2[!sexcheck2$sex.x==sexcheck2$sex.y,]$sex.x)
View(sexcheck2[!sexcheck2$sex.x==sexcheck2$sex.y&sexcheck2$sex.x!=0,])
table(sexcheck2[!sexcheck2$sex.x==sexcheck2$sex.y&sexcheck2$sex.x!=0,]$sex.x)
to_remove_sex_check<-as.data.frame(cbind(sexcheck2[!sexcheck2$sex.x==sexcheck2$sex.y,]$family.ID, sexcheck2[!sexcheck2$sex.x==sexcheck2$sex.y,]$sample.ID))
##14 outliers
colnames(to_remove_sex_check)<-c("FID","IID")
to_remove_sex_check<-to_remove_sex_check[!is.na(to_remove_sex_check$FID),]
write.table(to_remove_sex_check,"/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/Believe_geno_subset_to_remove_sex_check_plinkOK_HWE15.txt",sep='\t', row.names=FALSE,
            #col.names=FALSE,
            quote=FALSE, eol="\n")


#################################################
####Computation of PCs###########################
#################################################
library(bigsnpr)
# path_to_plink<-"/home/solene.cadiou/plink2.sh"
nc<-nb_cores()
# nc<-1
#######from genotype subset
##for the first call
obj.bed<-snp_readBed("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_correcthet.bed")
##after first time, use instead snp_attach
obj.bed<-snp_attach("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_correcthet.rds")
#####


##filter out some variants that are highly associated with population structure, e.g. as performed in the UK Biobank (Bycroft et al., 2018).##
###this step should be performed before relatedness computation
obj.bigsnp<-snp_attach("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/Solene_Believe_test/BELIEVE_genotype_forPCs_123456_HWE15_correcthet.rds")
G <- obj.bigsnp$genotypes
CHR <- obj.bigsnp$map$CHR
POS <- obj.bigsnp$map$POS

ind_chip <- 1:nrow(obj.bigsnp$map)
obj.bed <- bed("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_correcthet.bed")
##first round of PCS to identify vvariants
obj.svd<-bed_autoSVD(obj.bed, k = 5,
                     max.iter = 1)
plot(obj.svd)
plot(obj.svd, type = "scores")
plot(obj.svd, type = "scores", scores = 1:5)
plot(obj.svd, type = "scores", scores = 4:5)

# Get variants associated with pop struct, discard them and write bed file
# debugonce(snp_pcadapt)
# obj.pcadapt2 <- snp_pcadapt(G, obj.svd$u, ncores = nc)
# pcadapt0(G, obj.svd$u, ind.row, ind.col, ncores)
obj.pcadapt <- bed_pcadapt(obj.bed, obj.svd$u, ncores = nc)
plot(obj.pcadapt, type = "Manhattan")
plot(obj.pcadapt, type = "Q-Q")

length(ind_keep <- ind_chip[predict(obj.pcadapt, log10 = FALSE) > 0.05])
obj.bigsnp2 <- obj.bigsnp
##save the bed object without the variants associated to Pop
snp_writeBed(obj.bigsnp2, bedfile = "/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_correcthet_pcaaapt.bed", 
             ind.col = ind_keep)

# Compute the relatedness without these variants

rel <- runonce::save_run(
  snp_plinkKINGQC(path_to_plink, bedfile.in = "/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_correcthet_pcaaapt.bed",
                  thr.king = 2^-4.5, make.bed = FALSE, ncores =nc),
  file = "/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_correcthet_pcaaapt_rel.rds"
)
dim(rel)  
hist(rel$KINSHIP, breaks = 100); abline(v = 2^-(1.5:4.5), col = "red")
hist(log2(rel$KINSHIP), "FD"); abline(v = c(-4.5, -3.5), col = "red")

#########
#Remove related individuals
ind.rel <- match(c(rel$IID1, rel$IID2), obj.bed$fam$sample.ID)

head(ind.rel)
ind.norel <- rows_along(obj.bed)[-ind.rel]

length(ind.norel) ##6813
head(ind.norel)
##################
###Computation of PCs without these individuals
vobj.svd <- runonce::save_run(
  bed_autoSVD(obj.bed, k = 20,ind.row=ind.norel),
  file = "/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_correcthet_norel_pcadapt_PCA_removing_all.rds")
plot(vobj.svd)
plot(vobj.svd, type = "scores", scores = 1:20, coeff = 0.5)
plot(vobj.svd, type = "loadings", loadings = 1:20, coeff = 0.5)
PC <- predict(vobj.svd)

# saveRDS(PC,"PC_subset_genotype_123456789")
ldist <- log(bigutilsr::dist_ogk(PC))
# outlierness_on_subset<-cbind(IID=id,dist=ldist)
# outlierness_on_subset<-as.data.frame(outlierness_on_subset)
# outlier_id<-which(outlierness_on_subset$dist>=7.7)
hist(ldist, "FD")
plot_grid2(plotlist = lapply(1:4, function(k) {
  k1 <- 2 * k - 1
  k2 <- 2 * k
  qplot(PC[, k1], PC[, k2], color = ldist, size = I(2)) +
    scale_color_viridis_c() +
    theme_bigstatsr(0.6) +
    labs(x = paste0("PC", k1), y = paste0("PC", k2), color = "log-distance") +
    coord_equal()
}), nrow = 2, legend_ratio = 0.2, title_ratio = 0)

#Project PCs to remaining individuals 
PCs <- matrix(NA, nrow(obj.bed), ncol(vobj.svd$u))
dim(PCs)
head(PCs)

PCs[ind.norel, ] <- predict(vobj.svd)
head(PCs)

proj <- bed_projectSelfPCA(vobj.svd, obj.bed,ind.row = rows_along(obj.bed)[-ind.norel], ncores = 1)

PCs[-ind.norel, ] <- proj$OADP_proj
head(PCs)

sum(is.na(PCs))

#The order of the PCs should follow that of the fam file
colnames_list<-paste0('PC',1:20)
colnames(PCs)<-c(colnames_list)

PCs<-as.data.frame(PCs)
PCs$FID<-obj.bed$fam$family.ID
PCs$IID<-obj.bed$fam$sample.ID
PCs<- PCs %>% select(FID, IID, everything())
rel_TF<-rep(2,nrow(PCs))
rel_TF[ind.norel]<-rep(1,length(ind.norel))
plot(PCs$PC1,PCs$PC2, col=rel_TF)
plot(PCs$PC3,PCs$PC4, col=rel_TF)
plot(PCs$PC5,PCs$PC6, col=rel_TF)
plot(PCs$PC7,PCs$PC8, col=rel_TF)
plot(PCs$PC9,PCs$PC10, col=rel_TF)
plot(PCs$PC11,PCs$PC12, col=rel_TF)
##at this step we can decide to remove :
##1. the 2 outliers points in PC2 (identified as South Asian, see below)
##2. the clusters of PC3 and 4
##3. the clusters of PC5 and 6

##1. 2 outliers in PC2
outliers_PC2<-PCs[PCs$PC2>80,c("FID","IID")]
length(outliers_PC2$IID)
##2. outliers PC3
outliers_PC3<-PCs[PCs$PC3>80,c("FID","IID")]
length(outliers_PC3$IID)
##3. outliers PC5/PC6
outliers_PC5_6<-PCs[PCs$PC5>(30)&PCs$PC6<(-20),c("FID","IID")]
length(outliers_PC5_6$IID)
outliers<-c(outliers_PC2$IID,outliers_PC3$IID,outliers_PC5_6$IID)

PCs_WO<-PCs[!PCs$IID%in%outliers,]
plot(PCs_WO$PC1,PCs_WO$PC2)
plot(PCs_WO$PC3,PCs_WO$PC4)
plot(PCs_WO$PC5,PCs_WO$PC6)
plot(PCs_WO$PC7,PCs_WO$PC8)
plot(PCs_WO$PC9,PCs_WO$PC10)
plot(PCs_WO$PC11,PCs_WO$PC12)

##save list of samples to keep and outliers
ind.outliers<- which(obj.bed$fam$sample.ID%in%outliers)
ind.noout <- rows_along(obj.bed)[-ind.outliers]
# id_final<-obj.bed$fam$sample.ID[ind.noout]
# id_final_PCs_relall<-as.data.frame(cbind(FID=id_final,IID=id_final))
id_outliers_PCs_relall<-rbind(outliers_PC2,outliers_PC3,outliers_PC5_6)
table(duplicated(id_outliers_PCs_relall))
##no duplicate
write.table(id_outliers_PCs_relall,file="/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/GT_BELIEVE_subset_outliers_afterPCS_relall_correcthet.txt", sep='\t', row.names=FALSE,
            #col.names=FALSE,
            quote=FALSE, eol="\n")
id_final_PCs_relall<-PCs[!PCs$IID%in%id_outliers_PCs_relall$IID,c("FID","IID")]
dim(id_final_PCs_relall)
write.table(id_final_PCs_relall,file="/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/GT_BELIEVE_subset_to_keep_afterPCS_relall_correcthet.txt", sep='\t', row.names=FALSE,
            #col.names=FALSE,
            quote=FALSE, eol="\n")


#####################################################
######recomputation of PCs without these individuals
obj.bigsnp<-snp_attach("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/Solene_Believe_test/BELIEVE_genotype_forPCs_123456_HWE15_correcthet.rds")
G <- obj.bigsnp$genotypes
CHR <- obj.bigsnp$map$CHR
POS <- obj.bigsnp$map$POS

ind_chip <- 1:nrow(obj.bigsnp$map)
obj.bed <- bed("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_correcthet.bed")
obj.svd2<-bed_autoSVD(obj.bed, k = 5,
                      max.iter = 1, ind.row = ind.noout)
plot(obj.svd2)
plot(obj.svd2, type = "scores")
plot(obj.svd2, type = "scores", scores = 1:5)
plot(obj.svd2, type = "scores", scores = 4:5)

# Get variants associated with pop struct, discard them and write bed file
obj.pcadapt2 <- bed_pcadapt(obj.bed, obj.svd2$u, ind.row=ind.noout,ncores = 1)
plot(obj.pcadapt2, type = "Manhattan")
plot(obj.pcadapt2, type = "Q-Q")
length(ind_keep <- ind_chip[predict(obj.pcadapt2, log10 = FALSE) > 0.05])
obj.bigsnp3 <- obj.bigsnp
# obj.bigsnp2$map <- dplyr::transmute(
#   obj.bigsnp$map, chromosome = CHR, marker.ID = SNP, genetic.dist = 0,
#   physical.pos = POS, allele1 = a1, allele2 = a2)
snp_writeBed(obj.bigsnp3, bedfile = "/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_pcaaapt_WO_removing_all.bed", ind.col = ind_keep)

# Compute the relatedness

rel <- runonce::save_run(
  snp_plinkKINGQC(path_to_plink, bedfile.in = "/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/pQTL/Solene_Believe_test/BELIEVE_genotype_forPCs_123456_HWE15_pcaaapt_WO_removing_all.bed",
                  thr.king = 2^-4.5, make.bed = FALSE, ncores =nc),
  file = "/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_pcaaapt_rel_WO_removing_all.rds"
)
dim(rel)  
hist(rel$KINSHIP, breaks = 100); abline(v = 2^-(1.5:4.5), col = "red")
hist(log2(rel$KINSHIP), "FD"); abline(v = c(-4.5, -3.5), col = "red")


rel2 <- rbind(data.frame(IID = rel$IID1, K = rel$KINSHIP),
              data.frame(IID = rel$IID2, K = rel$KINSHIP)) %>%
  group_by(IID) %>%
  summarise(sum_K = sum(K))

rel3 <- subset(rel2, sum_K > 2^-3.5)
hist(rel2$sum_K, breaks = 100)
is_rel <- obj.bigsnp$fam$sample.ID %in% rel3$IID
sum(!is_rel)  # 7171
ind_to_keep<-which(!is_rel)
ind_to_keep2<-intersect(ind_to_keep,ind.noout)
obj.bed <- bed("/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_correcthet.bed")
id<-obj.bed$fam$sample.ID[ind_to_keep]
genid_subset_norel<-substr(id,19,27)

vobj.svd2 <- runonce::save_run(
  bed_autoSVD(obj.bed, k = 20,ind.row=ind.norel2),
  file = "/exchange/healthds/pQTL/BELIEVE/Genetic_QC_files/Intermediate_files_for_genetic_QC/BELIEVE_genotype_forPCs_123456_HWE15_norel_pcadapt_PCA_WO_removing_all.rds")
plot(vobj.svd2)
plot(vobj.svd2, type = "scores", scores = 1:20, coeff = 0.5)
plot(vobj.svd2, type = "loadings", loadings = 1:20, coeff = 0.5)
