###Import Libraries
remotes::install_github("privefl/bigstatsr")
remotes::install_github("privefl/bigsnpr")
install.packages("tidyverse")

library(bigsnpr)
library(bigstatsr)
library(bigreadr)
library(tidyverse)

library(ggplot2)
library(haven)
library(dplyr)
library(FactoMineR)
library(R.utils)
require(pgenlibr)
require(data.table)

devtools::source_url(
  "https://raw.githubusercontent.com/privefl/paper4-bedpca/master/code/plot_grid2.R")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("snpStats")



######## A.1. Compute common ID with proteomics dataset
## import samples present in the proteomics
path_prot <- "/center/healthds/pQTL/INTERVAL/cleaned_INTERVAL.Rds"
prot_data <- readRDS(path_prot, refhook = NULL)
prot_data <- data.frame(prot_data$imputed_cleaned_dataset) #9657
prot_ids <- prot_data[prot_data$SampleType == 'Sample', ]$SampleId

## import genomic non imputed data
fam_file<-read.table("/processing_data/shared_datasets/plasma_proteome/interval/genotypes/interval_qced_24.8.18.fam",sep="",header=F)
gen_ids <- as.factor(fam_file$V1) #42396
prot_ids <- as.factor(prot_ids)

id_conversion <- read.csv("/processing_data/shared_datasets/plasma_proteome/interval/phenotypes/INTERVAL_OmicsMap_20221221.csv", )
id_conversion <- id_conversion[,c(1,2,4,5,16)] 
id_conversion <- id_conversion %>% drop_na(Soma7000_RAW) #9769
id_conversion <- id_conversion[,c(1,3,5)] #9769
id_conversion <- id_conversion %>% drop_na() #9443
# the rest are NA in the conversion file

gen_ids_conv <- as.factor(id_conversion$Affymetrix_gwasQC_bl)
common_gen_ids <- Reduce(intersect, list(levels(gen_ids),levels(gen_ids_conv))) #9443
# we filter based on the proteomics data that we have

prot_ids_conv <- as.factor(id_conversion$Soma7000_RAW)
prot_ids_conv_str <- lapply(levels(prot_ids_conv), function(x) {if(nchar(x)<10) x= paste("0",x,sep="") else x})
common_prot_ids <- Reduce(intersect, list(levels(prot_ids),prot_ids_conv_str)) #9327

common_ids <- id_conversion[id_conversion$Soma7000_RAW %in% as.numeric(common_prot_ids),] #9327

# convert to string
common_ids$Affymetrix_gwasQC_bl <- as.character(common_ids$Affymetrix_gwasQC_bl)


# set samples to exclude
extra_ids = c('110017048382', '110016204883' ,'110008981087')

# exclude samples
common_ids = common_ids %>% filter(!Affymetrix_gwasQC_bl %in% extra_ids)
common_ids$Affymetrix_gwasQC_bl <- as.numeric(common_ids$Affymetrix_gwasQC_bl)

write.csv(common_ids, "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/common_ids_gen_prot.csv", row.names=FALSE)
# 9327 common ids between genomic and proteomics data


uncommon_ids <- id_conversion[!(id_conversion$Soma7000_RAW %in% as.numeric(common_prot_ids)),]
uncommon_ids$Affymetrix_gwasQC_bl = as.character(uncommon_ids$Affymetrix_gwasQC_bl)
uncommon_ids_check = common_ids %>% filter(Affymetrix_gwasQC_bl %in% extra_ids)
uncommon_ids$Affymetrix_gwasQC_bl = as.numeric(uncommon_ids$Affymetrix_gwasQC_bl)
# 116 uncommon ids

write.csv(uncommon_ids, "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/non_common_ids_gen_prot.csv", row.names=FALSE)

write.table(common_ids$Affymetrix_gwasQC_bl, "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step1/Samples to remove/common_ID.txt", sep = "\n", row.names = FALSE)

######## A.1. Extract common ID in genotype files in PLINK
#SRC_DIR=/processing_data/shared_datasets/plasma_proteome/interval/genotypes
#OUT_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step1/Common_ID/
#plink2 --bfile $SRC_DIR/interval_qced_24.8.18 --keep-fam $OUT_DIR/common_ID.txt --make-bed --out $OUT_DIR/interval_qced_24.8.18_restricted
#plink2 --bfile $SRC_DIR/merged_imputation --keep-fam $OUT_DIR/common_ID.txt --make-bed --out $OUT_DIR/merged_imputation_restricted

######## A.1. Extract common ID in imputed files in PLINK
#SRC_DIR=/processing_data/shared_datasets/plasma_proteome/interval/genotypes/imputed/pgen
#OUT_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step1/Common_ID/
#list=$(ls $SRC_DIR/*.pgen)
#for i in $(seq 2 22); do
# plink2 --pfile $SRC_DIR/impute_dedup_${i}_interval --keep-fam $OUT_DIR/common_ID.txt --make-pgen --out $OUT_DIR/pgen/chr${i}
#done

######## A.3. Compute heterozigosity in PLINK
#SRC_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step1/Common_ID
#OUT_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step3
#plink2 --bfile $SRC_DIR/merged_imputation_restricted --missing --out $OUT_DIR/missing_info
#plink2 --bfile $SRC_DIR/merged_imputation_restricted --not-chr X Y XY --maf 0.01 --geno 0.1 --mind 0.1 --hwe 1e-15 --make-bed --out $OUT_DIR/merged_imputation_restricted_QC
#plink2 --bfile $OUT_DIR/merged_imputation_restricted_QC --het --out $OUT_DIR/merged_imputation_heted_imputation_het

######## A.3. Identify heterozigosity outliers
het<-data.table::fread("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step3/merged_imputation_heted_imputation_het.het",
                       sep = "\t",header = TRUE, verbose = TRUE)
View(het)
m <- mean(het$F) # Calculate the mean  
s <- sd(het$F) # Calculate the SD
valid <- (het$F <= m+3*s & het$F >= m-3*s) 
# Get any samples with F coefficient within 3 SD of the population mean
table(valid)
##excluse 72 individuals
hist(het$F,100)
abline(v=m+3*s)
abline(v=m-3*s)
to_remove_het<-as.data.frame(cbind(FID=het$`#FID`[!valid],IID=het$IID[!valid]))
write.table(to_remove_het,file="/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step3/Interval_subset_to_remove_het_HWE15.txt", sep='\t', row.names=FALSE,
            col.names=FALSE, quote=FALSE, eol="\n")
##72


######## A.3. Remove heterozigosity outliers in PLINK and Create the bed file for PC computation
#SRC_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step1/Common_ID
#OUT_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step3
#plink2 --bfile $SRC_DIR/merged_imputation_restricted --remove $OUT_DIR/Interval_subset_to_remove_het_HWE15.txt --make-bed --out $OUT_DIR/merged_imputation_restricted_no_het_out
#plink2 --bfile $OUT_DIR/merged_imputation_restricted_no_het_out --not-chr X Y XY --maf 0.01 --geno 0.1 --mind 0.1 --hwe 1e-15 --make-bed --out $OUT_DIR/merged_imputation_PC

######## A.5.	Identification of samples with high relatedness to remove them in PCA
path_to_plink<-"/home/alessia.mapelli/plink2.sh"
nc<-nb_cores()
# nc<-1

#######from genotype subset
##for the first call
obj.bed<-snp_readBed("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step3/merged_imputation_PC.bed")
##after first time, use instead snp_attach
obj.bigsnp<-snp_attach("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step3/merged_imputation_PC.rds")
#####


##filter out some variants that are highly associated with population structure, e.g. as performed in the UK Biobank (Bycroft et al., 2018).##
###this step should be performed before relatedness computation
G <- obj.bigsnp$genotypes
CHR <- obj.bigsnp$map$CHR
POS <- obj.bigsnp$map$POS

ind_chip <- 1:nrow(obj.bigsnp$map)
obj.bed <- bed("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step3/merged_imputation_PC.bed")

##first round of PCS to identify variants
obj.svd<-bed_autoSVD(obj.bed, k = 5,
                     max.iter = 1)
# Fast truncated SVD with initial pruning and that iteratively removes long-range LD regions.
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
#Method to detect genetic markers involved in biological adaptation. 
#This provides a statistical tool for outlier detection based on Principal Component Analysis. 
#This corresponds to the statistic based on mahalanobis distance, as implemented in package pcadapt.

length(ind_keep <- ind_chip[predict(obj.pcadapt, log10 = FALSE) > 0.05])
obj.bigsnp2 <- obj.bigsnp
##save the bed object without the variants associated to Pop
snp_writeBed(obj.bigsnp2, bedfile = "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step5/merged_imputation_PC_pcaaapt.bed", 
             ind.col = ind_keep)
obj.bigsnp$fam$sample.ID <- as.character(obj.bigsnp$fam$sample.ID)

# Compute the relatedness without these variants
rel <- runonce::save_run(
  snp_plinkKINGQC(path_to_plink, bedfile.in = "//group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step5/merged_imputation_PC_pcaaapt.bed",
                  thr.king = 2^-4.5, make.bed = FALSE, ncores =nc),
  file = "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step5/merged_imputation_PC_pcaaapt_rel.rds"
)
dim(rel)  
hist(rel$KINSHIP, breaks = 100); abline(v = 2^-(1.5:4.5), col = "red")
hist(log2(rel$KINSHIP), "FD"); abline(v = c(-4.5, -3.5), col = "red")


##optimization of the list of samples to remove
library(dplyr)
rel2 <- rbind(data.frame(IID = rel$IID1, K = rel$KINSHIP),
              data.frame(IID = rel$IID2, K = rel$KINSHIP)) %>%
  group_by(IID) %>%
  summarise(sum_K = sum(K))

rel3 <- subset(rel2, sum_K > 2^-3.5)
hist(rel2$sum_K, breaks = 100)
sample_ids <- as.character(obj.bigsnp$fam$sample.ID)
rel_ids <- as.character(as.numeric(rel3$IID))
is_rel <- sample_ids %in% rel_ids
sum(!is_rel)  # 9196
ind_to_keep<-which(!is_rel)
#obj.bed <- bed("/group/diangelantonio/users/Solene/pQTL/Solene_Believe_test/BELIEVE_genotype_forPCs_123456_HWE15.bed")
# # obj.bed2 <- bed("/center/healthds/pQTL/Solene_Believe_test/BELIEVE_genotype_1234567_norel.bed")
# id<-obj.bed$fam$sample.ID[ind_to_keep]
# genid_subset_norel<-substr(id,19,27)

#########
#Remove all pairs of related individuals
ind.rel <- match(c(rel$IID1, rel$IID2), obj.bed$fam$sample.ID)

head(ind.rel)
length(ind.rel) #120, unique 117

ind.norel <- rows_along(obj.bed)[-ind.rel]
length(ind.norel) #9138
head(ind.norel)

snp_writeBed(obj.bigsnp2, bedfile = "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step5/merged_imputation_PC_pcaaapt_norel_r.bed", 
             ind.row = ind.norel, ind.col = ind_keep)


## Check with Plink2 that does the same and remove the 59 individuals
# obj.bed <- bed("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step5/merged_imputation_PC_pcaaapt_rel.bed")
# 9196 individuals

######## A.6.	Outliers for PCA

###Computation of PCs without related individuals
vobj.svd <- runonce::save_run(
  bed_autoSVD(obj.bed, k = 20, ind.row=ind.norel),
  file = "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step5/merged_imputation_PC_pcaaapt_removing_all.rds")
dev.off()
plot(vobj.svd)
plot(vobj.svd, type = "scores", scores = 1:8, coeff = 0.5)
plot(vobj.svd, type = "loadings", loadings = 1:8, coeff = 0.5)
PC <- predict(vobj.svd)

ldist <- log(bigutilsr::dist_ogk(PC))
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

#The order of PCs follows the fam file
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
PCs$rel <- rel_TF
save(PCs, file = "PCs.Rda")

load(file='PCs.Rda')

p<-ggplot(PCs,aes(x=PC1,y=PC2, colour=rel_TF))+
  geom_point()
p

# ----------------------------------------------------------
# COMPARE TO THOSE COMPUTED ON THE WHOLE POP
pca_all_int <- read.table("/processing_data/shared_datasets/plasma_proteome/interval/genotypes/annot_INT_50PCs_pcs.txt", header = TRUE)
common_ids <- read.csv("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step1/Common_ID/common_ids_gen_prot.csv")

present_in_prot <- common_ids$Affymetrix_gwasQC_bl
pca_all_int$prot_pres <- 0
pca_all_int[pca_all_int$ID %in% present_in_prot, ]$prot_pres<- 1

p<-ggplot(pca_all_int,aes(x=PC_1,y=PC_2, colour=prot_pres))+
  geom_point()
p
p<-ggplot(pca_all_int,aes(x=PC_3,y=PC_4, colour=prot_pres))+
  geom_point()
p



############ Impute ancestry #################################################
all_freq <- bigreadr::fread2(
  runonce::download_file(
    "https://figshare.com/ndownloader/files/38019027",  # for the tutorial (46 MB)
    # "https://figshare.com/ndownloader/files/31620968",  # for real analyses (849 MB)
    dir = "tmp-data", fname = "ref_freqs.csv.gz"))

projection <- bigreadr::fread2(
  runonce::download_file(
    "https://figshare.com/ndownloader/files/38019024",  # for the tutorial (44 MB)
    # "https://figshare.com/ndownloader/files/31620953",  # for real analyses (847 MB)
    dir = "tmp-data", fname = "projection.csv.gz"))

# coefficients to correct for overfitting of PCA
correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)

library(dplyr)
matched <- obj.bed$map %>%
  transmute(chr = chromosome, pos = physical.pos, a1 = allele1, a0 = allele2) %>% 
  mutate(beta = 1) %>%
  snp_match(all_freq[1:5], match.min.prop= 0.11) %>%
  print()
# <0.15 matched

# further subsetting on missing values
counts <- bed_counts(obj.bed, ind.col = matched$`_NUM_ID_.ss`)
hist(nb_na <- counts[4, ])
ind <- which(counts[4, ] < (nrow(obj.bed) * 0.05))  
length(ind)

# project individuals (divided by 2) onto the PC space
PROJ <- as.matrix(projection[matched$`_NUM_ID_`[ind], -(1:5)])
fun <- function(x) {
  bed_prodVec(obj.bed, x, ind.col = matched$`_NUM_ID_.ss`[ind],
              # scaling to get G if beta = 1 and (2 - G) if beta = -1
              center = 1 - matched$beta[ind],
              scale = matched$beta[ind]) }
all_proj <- apply(sweep(PROJ, 2, correction / 2, '*'), 2, fun)

all_centers <- crossprod(as.matrix(all_freq[matched$`_NUM_ID_`[ind], -(1:5)]), PROJ)
all_sq_dist <- apply(all_centers, 1, function(one_center) {
  rowSums(sweep(all_proj, 2, one_center, '-')^2)
})

THR <- 0.001  # you can adjust this threshold
thr_sq_dist <- max(dist(all_centers)^2) * THR / 0.16

group <- colnames(all_freq)[-(1:5)]
group[group %in% c("Scandinavia", "United Kingdom", "Ireland")]   <- "Europe (North West)"
group[group %in% c("Europe (South East)", "Europe (North East)")] <- "Europe (East)"

cluster <- apply(all_sq_dist, 1, function(sq_dist) {
  ind <- which.min(sq_dist)
  if (sq_dist[ind] < thr_sq_dist) group[ind] else NA
})

table(cluster, exclude = NULL)  # 209 NA

p<-ggplot(PCs,aes(x=PC1,y=PC2, colour=cluster))+
  geom_point()
p


p<-ggplot(PCs,aes(x=PC3,y=PC4, colour=cluster))+
  geom_point()
p

PCs$imputed_anchestry <- as.factor(cluster)
summary(PCs)
save(PCs, file = "PCs_imputed_anc.Rda")

################ Exploit the anchestry in INTERVAL data #####################
anchestry <- read.csv('/processing_data/shared_datasets/plasma_proteome/interval/phenotypes/INTERVALdata_21DEC2022.csv')
conversion <- read.csv('/processing_data/shared_datasets/plasma_proteome/interval/phenotypes/INTERVAL_OmicsMap_20221221.csv',)
#load(file = '/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step5/PCs_imputed_anc.Rda')
ethinic <- anchestry[,c(1,4)]
table(is.na(ethinic$ethnicPulse))
head(ethinic)
ethinic$ethnicPulse <- as.factor(ethinic$ethnicPulse)
true_anchestry_grouped <- lapply(anchestry$ethnicPulse, function(X){
  if(X == "" || X == "Not Disclosed"){X= "Unknown"} else {X}
})
ethinic$ethnicPulse_grouoped <- as.factor(unlist(true_anchestry_grouped))
summary(ethinic)
summary(conversion)
summary(PCs)
conversion <- conversion %>% drop_na(Affymetrix_gwasQC_bl)
gen_ids_conv <- as.factor(conversion$Affymetrix_gwasQC_bl)
gen_ids <- as.factor(PCs$FID)
common_gen_ids <- Reduce(intersect, list(levels(gen_ids),levels(gen_ids_conv))) #9255

common_ids <- conversion[conversion$Affymetrix_gwasQC_bl %in% as.numeric(common_gen_ids),] #9255
common_ids <- common_ids[c(1,4)]
df_merge <- merge(common_ids,ethinic,by="identifier")
df_merge$Affymetrix_gwasQC_bl <- as.factor(df_merge$Affymetrix_gwasQC_bl)
PCs$FID <- as.factor(PCs$FID)
summary(df_merge)
summary(PCs)
df_merge <- merge(x=df_merge,y=PCs,by.x="Affymetrix_gwasQC_bl", by.y= "FID")
col <- colnames(df_merge)
col[1] <- "FID"
col[3] <- "true_anchestry"
col[4] <- "true_anchestry_grouped"
colnames(df_merge) <- col
summary(df_merge)

write.csv(df_merge, "PCs_all_anc.csv")
save(df_merge, file = "PCs_all_anc.Rda")

df_merge <- read.csv("PCs_all_anc.csv")

p<-df_merge %>% filter(!true_anchestry_grouped == "Eng/W/Scot/NI/Brit") %>% ggplot(aes(x=PC3,y=PC4, colour=true_anchestry_grouped))+
  geom_point() + scale_color_brewer(palette="Paired")
p

p<-ggplot(df_merge,aes(x=PC3,y=PC4, colour=true_anchestry_grouped))+
  geom_point() + scale_color_brewer(palette="Paired")
p

table(df_merge$true_anchestry_grouped)

######## B.	Variants with low imputation quality

######## B.1. Recompute summary metrics within the sample with proteomics data in PLINK
#OUT_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepB/pgen_restricted_sample
#FAM_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step5-6
#SRC_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step1/Common_ID/pgen

#for i in $(seq 1 22); do
# plink2 --pfile $SRC_DIR/chr${i} --threads 16 --memory 16384 --keep-fam $FAM_DIR/merged_imputation_PC_pcaaapt_rel.fam --export bgen-1.2 --out $OUT_DIR/chr${i}
# qctool -g $OUT_DIR/chr${i}.bgen -snp-stats -osnp $OUT_DIR/snp-stats_chr${i}.txt
#done

######## B.2. o	Exclude variants with --mac 20 --hwe 1e-15 info_score > 0.7
path_to_snpstat_new <- "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepB/pgen_restricted_sample"
path_to_save <- "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepB/SummaryQC"
intital_n_var <- 0
sum_keeped <- 0
sum_keeped_over_thr <- 0
i <-1
for(i in seq(1,22)){
  chr <- paste("snp-stats_chr",i,".txt", sep="")
  chr_old <- read.table(paste(path_to_snpstat, "/impute_", i,"_interval.snpstats", sep= ""), header = T)
  intital_n_var <- intital_n_var + dim(chr_old)[1]
  info <- read.table(paste(path_to_snpstat_new, chr, sep= "/"), header = T)
  summary(info)
  #save_png <- paste("snp-stats_chr",i,".png", sep="")
  #png(file=paste(path_to_save, save_png, sep= "/"),
  #    width=600, height=350)
  #plot(density(info$info), type="l", main=paste("Density info score CHR ", i, sep=""))
  #dev.off()
  sum_keeped <- sum_keeped + dim(info)[1]
  #chr_keeped_under_0.7 <- info[info$info < 0.7, c(2,14,17)]
  #write.table(chr_keeped_under_0.7, file=paste(path_to_save, "/keeped_under_0.7_chr_", i, sep= ""), row.names = FALSE)
  sum_keeped_over_thr <- sum_keeped_over_thr + dim(info[info$info > 0.7, ])[1]
  #write.table(info[info$info > 0.7, ], file=paste(path_to_save, "/keeped_over_0.7_chr_", i, sep= ""), row.names = FALSE)
  #write.table(info, file=paste(path_to_save, "/keeped_snp_chr_", i, sep= ""), row.names = FALSE)
  #write.table(info$rsid, file=paste(path_to_save, "/snp_id_chr_", i, sep= ""), row.names = FALSE , col.names = F)
  write.table(info[info$info > 0.7, ]$rsid, file=paste(path_to_save, "/keeped_snp_over_0.7_chr_", i, sep= ""), row.names = FALSE, col.names = F,quote=F)
  
}
c(intital_n_var,sum_keeped,sum_keeped_over_thr)
write.table(c(intital_n_var,sum_keeped,sum_keeped_over_thr), file=paste(path_to_save, "/snps_count", sep= ""), row.names = FALSE, col.names = F)

######## C. Preparation of the final dataset

######## C.3. Genotype files in PLINK
# SRC_DIR=/processing_data/shared_datasets/plasma_proteome/interval/genotypes
# OUT_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC
# ID_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step1/Common_ID/
#   FAM_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step3
# 
# 
# plink2 \
# --bfile $SRC_DIR/merged_imputation \
# --keep-fam $FAM_DIR/merged_imputation_restricted_no_het_out.fam \
# --not-chr X Y XY \
# --geno 0.1 \
# --mind 0.1 \
# --mac 20 \
# --hwe 1e-15 \
# --make-bed \
# --out $OUT_DIR/cleaned_genotype_INTERVAL



######## C.1. Imputed pgen in PLINK
# SRC_DIR=/processing_data/shared_datasets/plasma_proteome/interval/genotypes/imputed/pgen
# OUT_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC
# ID_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepB/SummaryQC/
#   for i in $(seq 1 22); do
# plink2 \
# --pfile $SRC_DIR/impute_dedup_${i}_interval \
# --keep-fam $OUT_DIR/cleaned_genotype_INTERVAL.fam \
# --extract $ID_DIR/keeped_snp_over_0.7_chr_${i} \
# --not-chr X Y XY \
# --geno 0.1 \
# --mind 0.1 \
# --mac 20 \
# --hwe 1e-15 \
# --make-pgen \
# --out $OUT_DIR/pgen/cleaned_imputed_INTERVAL_chr_${i}
# done

# source /center/healthds/singularity_functions
# cd $HOME
# SRC_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC/pgen
# OUT_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC
# for i in $(seq 1 22); do
# plink2 --pfile $SRC_DIR/cleaned_imputed_INTERVAL_chr_${i} --recode vcf bgz id-paste=iid --out $OUT_DIR/vcf/chr$i
# done
# files=$(ls $OUT_DIR/vcf/*.vcf.gz)
# bcftools concat $files -Oz -o $OUT_DIR/vcf/allchromosomes.vcf.gz
# plink2 --vcf $OUT_DIR/vcf/allchromosomes.vcf.gz --make-pgen --out $OUT_DIR/pgen/allchromosomes_imputed
# plink2 --pfile $OUT_DIR/pgen/allchromosomes_imputed --pgen-info


######## C.2. Imputed bgen in PLINK
# OUT_DIR=/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC
# for i in $(seq 1 22); do
# plink2 \
# --pfile $OUT_DIR/pgen/cleaned_imputed_INTERVAL_chr_${i} \
# --export bgen-1.2 \
# --out $OUT_DIR/bgen/cleaned_imputed_INTERVAL_chr_${i}
# done
# plink2 \
# --pfile $OUT_DIR/pgen/allchromosomes_imputed \
# --export bgen-1.2 \
# --out $OUT_DIR/bgen/allchromosomes_imputed

###### D. Compute the first 20 PCs
setwd("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC")

library(bigsnpr)
path_to_plink<-"/home/alessia.mapelli/plink2.sh"
nc<-nb_cores()
# nc<-1

#######from genotype subset
##for the first call
obj.bed<-snp_readBed("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC/cleaned_genotype_INTERVAL.bed")
##after first time, use instead snp_attach
obj.bigsnp<-snp_attach("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC/cleaned_genotype_INTERVAL.rds")
#####


##filter out some variants that are highly associated with population structure, e.g. as performed in the UK Biobank (Bycroft et al., 2018).##
###this step should be performed before relatedness computation
G <- obj.bigsnp$genotypes
CHR <- obj.bigsnp$map$CHR
POS <- obj.bigsnp$map$POS

ind_chip <- 1:nrow(obj.bigsnp$map)
obj.bed <- bed("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC/cleaned_genotype_INTERVAL.bed")

##first round of PCS to identify variants
obj.svd<-bed_autoSVD(obj.bed, k = 5,
                     max.iter = 1)
# Fast truncated SVD with initial pruning and that iteratively removes long-range LD regions.
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
#Method to detect genetic markers involved in biological adaptation. 
#This provides a statistical tool for outlier detection based on Principal Component Analysis. 
#This corresponds to the statistic based on mahalanobis distance, as implemented in package pcadapt.

length(ind_keep <- ind_chip[predict(obj.pcadapt, log10 = FALSE) > 0.05])
obj.bigsnp2 <- obj.bigsnp
##save the bed object without the variants associated to Pop
snp_writeBed(obj.bigsnp2, bedfile = "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC/cleaned_genotype_INTERVA_pcaaapt.bed", 
             ind.col = ind_keep)
obj.bigsnp$fam$sample.ID <- as.character(obj.bigsnp$fam$sample.ID)

# Compute the relatedness without these variants
rel <- runonce::save_run(
  snp_plinkKINGQC(path_to_plink, bedfile.in = "//group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC/cleaned_genotype_INTERVA_pcaaapt.bed",
                  thr.king = 2^-4.5, make.bed = FALSE, ncores =nc),
  file = "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC/cleaned_genotype_INTERVA_pcaaapt_rel.rds"
)
dim(rel)  
hist(rel$KINSHIP, breaks = 100); abline(v = 2^-(1.5:4.5), col = "red")
hist(log2(rel$KINSHIP), "FD"); abline(v = c(-4.5, -3.5), col = "red")


##optimization of the list of samples to remove
library(dplyr)
rel2 <- rbind(data.frame(IID = rel$IID1, K = rel$KINSHIP),
              data.frame(IID = rel$IID2, K = rel$KINSHIP)) %>%
  group_by(IID) %>%
  summarise(sum_K = sum(K))

rel3 <- subset(rel2, sum_K > 2^-3.5)
hist(rel2$sum_K, breaks = 100)
sample_ids <- as.character(obj.bigsnp$fam$sample.ID)
rel_ids <- as.character(as.numeric(rel3$IID))
is_rel <- sample_ids %in% rel_ids
sum(!is_rel)  # 9196
ind_to_keep<-which(!is_rel)
#obj.bed <- bed("/group/diangelantonio/users/Solene/pQTL/Solene_Believe_test/BELIEVE_genotype_forPCs_123456_HWE15.bed")
# # obj.bed2 <- bed("/center/healthds/pQTL/Solene_Believe_test/BELIEVE_genotype_1234567_norel.bed")
# id<-obj.bed$fam$sample.ID[ind_to_keep]
# genid_subset_norel<-substr(id,19,27)

#########
#Remove all pairs of related individuals
ind.rel <- match(c(rel$IID1, rel$IID2), obj.bed$fam$sample.ID)

head(ind.rel)
length(ind.rel) #120, unique 117

ind.norel <- rows_along(obj.bed)[-ind.rel]
length(ind.norel) #9138
head(ind.norel)

snp_writeBed(obj.bigsnp2, bedfile = "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC/cleaned_genotype_INTERVA_pcaaapt_norel_r.bed", 
             ind.row = ind.norel, ind.col = ind_keep)


## Check with Plink2 that does the same and remove the 59 individuals
# obj.bed <- bed("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step5/merged_imputation_PC_pcaaapt_rel.bed")
# 9196 individuals

##################
###Computation of PCs without these individuals
vobj.svd <- runonce::save_run(
  bed_autoSVD(obj.bed, k = 20, ind.row=ind.norel),
  file = "/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/StepC/cleaned_genotype_INTERVA_pcaaapt_removing_all.rds")
dev.off()
plot(vobj.svd)
plot(vobj.svd, type = "scores", scores = 1:8, coeff = 0.5)
plot(vobj.svd, type = "loadings", loadings = 1:8, coeff = 0.5)
PC <- predict(vobj.svd)

ldist <- log(bigutilsr::dist_ogk(PC))
hist(ldist, "FD")
plot_grid(plotlist = lapply(1:2, function(k) {
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

#The order of PCs follows the fam file
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
PCs$rel <- rel_TF
save(PCs, file = "PCs.Rda")

load(file='PCs.Rda')

dev.off()
p<-ggplot(PCs,aes(x=PC1,y=PC2, colour=rel))+
  geom_point()
p

write.csv(PCs, file = "PCs_rel.csv", row.names=FALSE)

# ----------------------------------------------------------
# COMPARE TO THOSE COMPUTED ON THE WHOLE POP
pca_all_int <- read.table("/processing_data/shared_datasets/plasma_proteome/interval/genotypes/annot_INT_50PCs_pcs.txt", header = TRUE)
common_ids <- read.csv("/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step1/Common_ID/common_ids_gen_prot.csv")

present_in_prot <- common_ids$Affymetrix_gwasQC_bl
pca_all_int$prot_pres <- 0
pca_all_int[pca_all_int$ID %in% present_in_prot, ]$prot_pres<- 1

p<-ggplot(pca_all_int,aes(x=PC_1,y=PC_2, colour=prot_pres))+
  geom_point()
p
p<-ggplot(pca_all_int,aes(x=PC_3,y=PC_4, colour=prot_pres))+
  geom_point()
p



############ Impute ancestry #################################################
all_freq <- bigreadr::fread2(
  runonce::download_file(
    "https://figshare.com/ndownloader/files/38019027",  # for the tutorial (46 MB)
    # "https://figshare.com/ndownloader/files/31620968",  # for real analyses (849 MB)
    dir = "tmp-data", fname = "ref_freqs.csv.gz"))

projection <- bigreadr::fread2(
  runonce::download_file(
    "https://figshare.com/ndownloader/files/38019024",  # for the tutorial (44 MB)
    # "https://figshare.com/ndownloader/files/31620953",  # for real analyses (847 MB)
    dir = "tmp-data", fname = "projection.csv.gz"))

# coefficients to correct for overfitting of PCA
correction <- c(1, 1, 1, 1.008, 1.021, 1.034, 1.052, 1.074, 1.099,
                1.123, 1.15, 1.195, 1.256, 1.321, 1.382, 1.443)

library(dplyr)
matched <- obj.bed$map %>%
  transmute(chr = chromosome, pos = physical.pos, a1 = allele1, a0 = allele2) %>% 
  mutate(beta = 1) %>%
  snp_match(all_freq[1:5], match.min.prop= 0.11) %>%
  print()
# <0.15 matched

# further subsetting on missing values
counts <- bed_counts(obj.bed, ind.col = matched$`_NUM_ID_.ss`)
hist(nb_na <- counts[4, ])
ind <- which(counts[4, ] < (nrow(obj.bed) * 0.05))  
length(ind)

# project individuals (divided by 2) onto the PC space
PROJ <- as.matrix(projection[matched$`_NUM_ID_`[ind], -(1:5)])
fun <- function(x) {
  bed_prodVec(obj.bed, x, ind.col = matched$`_NUM_ID_.ss`[ind],
              # scaling to get G if beta = 1 and (2 - G) if beta = -1
              center = 1 - matched$beta[ind],
              scale = matched$beta[ind]) }
all_proj <- apply(sweep(PROJ, 2, correction / 2, '*'), 2, fun)

all_centers <- crossprod(as.matrix(all_freq[matched$`_NUM_ID_`[ind], -(1:5)]), PROJ)
all_sq_dist <- apply(all_centers, 1, function(one_center) {
  rowSums(sweep(all_proj, 2, one_center, '-')^2)
})

THR <- 0.001  # you can adjust this threshold
thr_sq_dist <- max(dist(all_centers)^2) * THR / 0.16

group <- colnames(all_freq)[-(1:5)]
group[group %in% c("Scandinavia", "United Kingdom", "Ireland")]   <- "Europe (North West)"
group[group %in% c("Europe (South East)", "Europe (North East)")] <- "Europe (East)"

cluster <- apply(all_sq_dist, 1, function(sq_dist) {
  ind <- which.min(sq_dist)
  if (sq_dist[ind] < thr_sq_dist) group[ind] else NA
})

table(cluster, exclude = NULL)  # 209 NA

p<-ggplot(PCs,aes(x=PC1,y=PC2, colour=cluster))+
  geom_point()
p


p<-ggplot(PCs,aes(x=PC3,y=PC4, colour=cluster))+
  geom_point()
p

PCs$imputed_anchestry <- as.factor(cluster)
summary(PCs)
save(PCs, file = "PCs_imputed_anc.Rda")
write.csv(PCs, file = "PCs_imputed_anc.csv", row.names=FALSE)

################ Exploit the anchestry in INTERVAL data #####################
anchestry <- read.csv('/processing_data/shared_datasets/plasma_proteome/interval/phenotypes/INTERVALdata_21DEC2022.csv')
conversion <- read.csv('/processing_data/shared_datasets/plasma_proteome/interval/phenotypes/INTERVAL_OmicsMap_20221221.csv',)
#load(file = '/group/diangelantonio/users/alessia_mapelli/QC_gen_INTERVAL/QC_steps/Step5/PCs_imputed_anc.Rda')
ethinic <- anchestry[,c(1,4)]
table(is.na(ethinic$ethnicPulse))
head(ethinic)
ethinic$ethnicPulse <- as.factor(ethinic$ethnicPulse)
true_anchestry_grouped <- lapply(anchestry$ethnicPulse, function(X){
  if(X == "" || X == "Not Disclosed"){X= "Unknown"} else {X}
})
ethinic$ethnicPulse_grouoped <- as.factor(unlist(true_anchestry_grouped))
summary(ethinic)
summary(conversion)
summary(PCs)
conversion <- conversion %>% drop_na(Affymetrix_gwasQC_bl)
gen_ids_conv <- as.factor(conversion$Affymetrix_gwasQC_bl)
gen_ids <- as.factor(PCs$FID)
common_gen_ids <- Reduce(intersect, list(levels(gen_ids),levels(gen_ids_conv))) #9255

common_ids <- conversion[conversion$Affymetrix_gwasQC_bl %in% as.numeric(common_gen_ids),] #9255
common_ids <- common_ids[c(1,4)]
df_merge <- merge(common_ids,ethinic,by="identifier")
df_merge$Affymetrix_gwasQC_bl <- as.factor(df_merge$Affymetrix_gwasQC_bl)
PCs$FID <- as.factor(PCs$FID)
summary(df_merge)
summary(PCs)
df_merge <- merge(x=df_merge,y=PCs,by.x="Affymetrix_gwasQC_bl", by.y= "FID")
col <- colnames(df_merge)
col[1] <- "FID"
col[3] <- "true_anchestry"
col[4] <- "true_anchestry_grouped"
colnames(df_merge) <- col
summary(df_merge)

write.csv(df_merge, "PCs_all_anc.csv")
save(df_merge, file = "PCs_all_anc.Rda")

df_merge <- read.csv("PCs_all_anc.csv")

p<-df_merge %>% filter(!true_anchestry_grouped == "Eng/W/Scot/NI/Brit") %>% ggplot(aes(x=PC3,y=PC4, colour=true_anchestry_grouped))+
  geom_point() + scale_color_brewer(palette="Paired")
p

p<-ggplot(df_merge,aes(x=PC3,y=PC4, colour=true_anchestry_grouped))+
  geom_point() + scale_color_brewer(palette="Paired")
p

table(df_merge$true_anchestry_grouped)
