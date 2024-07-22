## /work/cgiambartolomei/SCRATCH/mvp_test/GTEXv8_COLOC_covid_cond.R
#args <- commandArgs(trailingOnly=TRUE)
#i = as.numeric(args[1]) # biom_file (4 files)
#j = as.numeric(args[2]) # tissues (49 tissues)
#ntissue = as.numeric(args[1])

wd = "/Users/c.giambartolomei/Library/CloudStorage/Dropbox/MVP/DANIELLE/"
setwd(wd)

file_paths = read.table("FILES/file_paths.txt", header=T)
genes_pos = read.table("/Users/c.giambartolomei/Library/CloudStorage/Dropbox/MVP/DANIELLE/FILES/gene_positions.txt", header=T, stringsAsFactors=F)
file_paths = cbind(file_paths, genes_pos)

file_paths$df2_file = gsub("/work/cgiambartolomei/SCRATCH/DANIELLE/","/Users/c.giambartolomei/Library/CloudStorage/Dropbox/MVP/DANIELLE/", file_paths$df2_file)
file_paths$df1_file = gsub("/work/cgiambartolomei/SCRATCH/DANIELLE/eQTL","/Users/c.giambartolomei/Library/CloudStorage/Dropbox/MVP/DANIELLE/lifted_from38_to37", file_paths$df1_file)
biom_file = file_paths$df2_file[1]


### OUTPUT
outfolder = paste(wd, "/out_conditional/", sep="")
if (!file.exists(outfolder)) dir.create(outfolder)

tempfolder = paste(outfolder, "temp_folder/", sep="")
if (!file.exists(tempfolder)) dir.create(tempfolder)

outfolder = paste(wd, "/GWAS_matched/", sep="")
if (!file.exists(outfolder)) dir.create(outfolder)

plotfolder = paste(wd, "/plots/", sep="")
if (!file.exists(plotfolder)) dir.create(plotfolder)

###########
### LIBRARIES
#############
library(data.table)
library(susieR)
#library(moloc)
library(coloc)
library(foreach)
library(doParallel)
library(GenomicRanges)
# library("Rmpfr")
library(Rfast)
# devtools::install_github('RfastOfficial/Rfast')
library(locuszoomr)

###########
### INPUT AND PATHS
#############
# only include SNPs above a MAF
maf_threshold = 0.001
p1_coloc = 1e-05
p2_coloc = 1e-05
p12_coloc = 1e-06

# With LD, must match by rsid since the LD info uses rsid

#path.plink <- "/work/cgiambartolomei/software/plink"
path.plink <- "/Applications/plink"
path.plink2 <- "/Applications/plink2"
# gcta_path = "/home2/cgiambartolomei/software/gcta_1.93.2beta/"
gcta_path = "/Users/c.giambartolomei/Library/CloudStorage/Dropbox/MVP/MVP_COVID_LIAM/hmgcr/gcta_1.93.2beta_mac/gcta64"

path.1kG <- "/work/gen_ref_datasets/plink_1KG_bypop/EUR/chr"
path.UKB <- "/home/cgiambartolomei/UKB_for_LD_EUR/sub/ukb.chr" # "ukb.chr20.impV3.HWE.EUR_sub.bim"
#path.GTEX <- "genotypes/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig_"
path.GTEX <- "genotypes/GTEx.WGS.838.passOnly.geno0.05.hwe0.00001.dbsnp.SNPsOnly.NoAmbig"




##############
### ANALYSIS
##############
#############
### eqtl data
#############

for (i in 1:nrow(file_paths)) {
        biom_file = file_paths$df2_file [ i ]
        eqtl_file = file_paths$df1_file [ i ]
        prefix = paste(basename(biom_file), basename(eqtl_file), sep="_")

        outfile = paste(outfolder, prefix, sep="")
        d = read.table(outfile, header=T,stringsAsFactors=F)
        ##########################
        # Extract LD, Freq, harmonize
        #############
        # LD info is recorded using rsID
        # add rsid to data from other dataset or fetch it
        chr = unique(d$CHR_eqtl)
        list.snps = d$rsID_eqtl
        list.reference.allele = data.frame(d$rsID_eqtl, d$A1_eqtl)

        # align both datasets to GWAS data alleles?
        # temp files
        prefix_temp = paste(prefix, "_temp", sep="")
        list.snps_file = paste(tempfolder, prefix_temp, "_list.snps", sep="")
        list.reference.allele_file = paste(tempfolder, prefix_temp, "_list.reference.allele", sep="")
        ld_file = paste(tempfolder, prefix_temp, "_ld", sep="")
        freq_file = paste(tempfolder, prefix_temp, "_freq", sep="")

        write.table(x=list.snps, file=list.snps_file, row.names = FALSE, col.names=FALSE, quote = FALSE, sep = '\t')
        write.table(x=list.reference.allele, file=list.reference.allele_file, row.names = FALSE, col.names=FALSE, quote = FALSE, sep = '\t')
        # could force Based on Christopher Chang's answer, to compute LD in Plink using a specific allele as a reference we can use:
        # plink --bfile 1000Genomes --r2 square --extract list.snps ---a2-allele list.reference.allele

        # Remove monomorphic SNPs in the dataset (those with a MAF = 0.0)
        cmd1 = paste(path.plink, " -bfile ", path.GTEX, " --maf ", maf_threshold, " --make-bed --r2 square --extract ", list.snps_file, " --a1-allele ", list.reference.allele_file, " --write-snplist --out ", ld_file, sep="")
#system(cmd1)

cmd1_error  <- tryCatch(system(cmd1), error=function(e) NULL )
if (cmd1_error!=0) {
        # could be because of duplicate IDs in the genotype data
        message("Plink command didn't work: try removing duplicate Ids?", file=stderr())
        dupl_file = paste(tempfolder, prefix_temp, "_dupl", sep="")

        #cmd_remove = paste(path.plink, " -bfile ", path.1kG, chr, " --list-duplicate-vars ids-only suppress-first --out ", dupl_file, sep="")
        #system(cmd_remove)
        # now try again
        #cmd1 = paste(path.plink, " -bfile ", path.1kG, chr, " --maf ", maf_threshold, " --make-bed --r2 square --extract ", list.snps_file, " --a1-allele ", list.reference.allele_file, " --write-snplist --out ", ld_file, " -exclude ", dupl_file, ".dupvar", sep="")
        #cmd1_error  <- tryCatch(system(cmd1), error=function(e) NULL )
        cmd_remove = paste(path.plink2, " -bfile ", path.GTEX, " --rm-dup ", " --maf ", maf_threshold, " --extract ", list.snps_file, " --write-snplist --out ", dupl_file, sep="")
        system(cmd_remove)

        if (file.exists(paste(dupl_file, ".rmdup.mismatch", sep=""))) {
        # remove duplicate snp from list of snps file
        dupl = scan(paste(dupl_file, ".rmdup.mismatch", sep=""),what="character")
        list.snps = list.snps[!(list.snps %in% dupl)]
        list.reference.allele = list.reference.allele[!(list.reference.allele[,1] %in% dupl),]
        write.table(x=list.snps, file=list.snps_file, row.names = FALSE, col.names=FALSE, quote = FALSE, sep = '\t')
        write.table(x=list.reference.allele, file=list.reference.allele_file, row.names = FALSE, col.names=FALSE, quote = FALSE, sep = '\t')

        # now try again
        cmd1 = paste(path.plink, " -bfile ", path.GTEX, " --maf ", maf_threshold, " --make-bed --r2 square --extract ", list.snps_file, " --a1-allele ", list.reference.allele_file, " --write-snplist --out ", ld_file, " -exclude ", dupl_file, ".rmdup.mismatch", sep="")
        cmd1_error  <- tryCatch(system(cmd1), error=function(e) NULL )
        }
}

if (!file.exists(paste(ld_file, ".ld", sep="")))  {
        print("###### !!!!! #######")
        print(prefix)
        message("The LD file for ", prefix,  " could not be created, possibly because Error: No variants remaining after main filters.")
cmd1 = paste(path.plink, " -bfile ", path.GTEX, " --maf ", maf_threshold, " --make-bed --r2 square --extract ", list.snps_file, " --a1-allele ", list.reference.allele_file, " --write-snplist --out ", ld_file, sep="")
cmd1_error  <- tryCatch(system(cmd1), error=function(e) NULL )
        print("###### !!!!! #######")
# remove all gcta files
do.call(file.remove, list(list.files(tempfolder, full.names = TRUE, pattern=prefix_temp)))
next()
}


cmd2=paste(path.plink," -bfile ", path.GTEX, " --maf ", maf_threshold, " --freq --extract ", list.snps_file, " --a1-allele ", list.reference.allele_file, " --out ", freq_file, sep="")
system(cmd2)

LD = as.matrix(fread(paste0(ld_file, ".ld"), data.table = F))
LD_snps = scan(paste0(ld_file, ".snplist"), what="character")
colnames(LD) <- rownames(LD) <- LD_snps

freq = fread(paste(freq_file, ".frq", sep=""), data.table = F)

x = d
names(x)[names(x) == 'rsID_eqtl'] <- 'SNPID'
x$Freq1_new = freq$MAF[match(x$SNPID, freq$SNP)]

# remove nan
ind_cols = colnames(LD)[which(is.nan(LD))]
ind_cols = unique(ind_cols)
LD = LD[,!(colnames(LD) %in% ind_cols)]
ind_rows = rownames(LD)[which(is.nan(LD))]
ind_rows = unique(ind_rows)
LD = LD[!(rownames(LD) %in% ind_rows),]
# ind=which(is.nan(LD))
toremove = unique(c(ind_cols, ind_rows))
x = x[!(x$SNPID %in% toremove),]

# make sure LD and data is in the same order
snps = colnames(LD)
x = x[match(snps, x$SNPID),]


snp_minp = d$rsID_eqtl[which.min(d$PVAL_eqtl)]
snp_minp_col = LD[,snp_minp]
d$LD_snp_minp = snp_minp_col[match(d$rsID_eqtl,names(snp_minp_col))]

# re-create data
#s_biom = unique(x$Ncases_biom/x$N_biom)
#dataset.biom = list(snp = x$SNPID, beta = x$BETA_biom, varbeta= (x$SE_biom)^2, s=s_biom, type = "cc", MAF=x$MAF_biom, N=x$N_biom, position = x$POS_biom, LD=LD)
#if (!is.null(check_dataset(dataset.biom,req="LD"))) stop("Problem with LD in biom data")

#dataset.eqtl = list(snp = x$SNPID, beta = x$BETA_eqtl, varbeta= (x$SE_eqtl)^2, type = "quant", MAF=x$MAF_eqtl, N=x$N_eqtl, position = x$POS_biom, LD = LD)
#if (!is.null(check_dataset(dataset.eqtl,req="LD"))) stop("Problem with LD in eqtl data")

#par(mfrow=c(2,1))
#plot_dataset(dataset.biom, main="GWAS")
#plot_dataset(dataset.eqtl, main="eQTL")

#dataset_biom_susie = runsusie(dataset.biom, R_ref=715, n=394846)
#dataset_eqtl_susie = runsusie(dataset.eqtl, n=715)

gap = 500000 # for visual comfort
library(EnsDb.Hsapiens.v75)
gwas = data.frame(chrom = d$CHR_biom, pos = d$POS_37_biom, rsid = d$rsID_eqtl, p = d$PVAL_biom, r2 = d$LD_snp_minp^2)
loc <- locus(data=gwas, ens_db = "EnsDb.Hsapiens.v75", seqname=17, xrange=c(min(gwas$pos)-gap, max(gwas$pos)+gap), LD = "r2")
#locus_plot(loc)

eqtl = data.frame(chrom = d$CHR_biom, pos = d$POS_37_biom, rsid = d$rsID_eqtl, p = d$PVAL_eqtl, r2 = d$LD_snp_minp^2)
loc2 <- locus(data=eqtl, ens_db = "EnsDb.Hsapiens.v75", seqname=17, xrange=c(min(gwas$pos)-gap, max(gwas$pos)+gap), LD = "r2")

plot_file = paste(plotfolder, prefix, ".png", sep="")
png(plot_file)
#pdf("myplot.pdf", width = 9, height = 6)
oldpar <- set_layers(2)
scatter_plot(loc, xticks = FALSE, legend_pos = 'topleft')
title("GWAS")
scatter_plot(loc2, legend_pos = NULL)
title("eQTL")
genetracks(loc, filter_gene_biotype = "protein_coding")
par(oldpar)  # revert par() settings

dev.off()

#par(mfrow=c(2,1))
#               locus_plot(loc, legend_pos = 'topleft', labels = "index", xticks = FALSE)
#               locus_plot(loc2, legend_pos = NULL,
#                          labels = "index", filter_gene_biotype = "protein_coding")
#dev.off()



plot_susie = FALSE
if (plot_susie) {
plot_file = paste(plotfolder, prefix, ".png", sep="")
png(plot_file)
#par(mfrow=c(2,1))
#plot_dataset(dataset.biom, susie_obj=dataset_biom_susie, main="GWAS")
#plot_dataset(dataset.eqtl, susie_obj=dataset_eqtl_susie, main="eQTL")
# sensitivity(susie_g,"H4 > 0.9",row=1,dataset1=dataset.biom,dataset2=dataset.eqtl)
# sensitivity(susie_g,"H4 > 0.9")
dev.off()

plot_error  <- tryCatch(sensitivity(susie_g,"H4 > 0.9"), error=function(e) NULL )
if (!is.null(plot_error)) {
plot2_file = paste(plotfolder, prefix, "COLOC.png", sep="")
png(plot2_file)
message("UHUUUUU FOUND A COLOC!!!!!")
sensitivity(susie_g,"H4 > 0.9")
dev.off()
}
}
