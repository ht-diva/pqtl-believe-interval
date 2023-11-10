#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("biomaRt")

library(biomaRt)
# super useful function 'separate_rows' from tidyr to separate rows by delimiter 
library(tidyr)

# Aptamer level annotations were created by mapping proteins to genomic coordinates using GENCODE (GRCh38), version 32 (Ensembl 98) [76].
# GENCODE: http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&c=chrX&g=knownGene

pr_file = "/center/healthds/pQTL/Reference_datasets_for_QC_proteomics/Somalogic/SomaScan_V4.1_7K_Annotated_Content_20210616.csv"
pr = read.csv(pr_file, header=T,stringsAsFactors=F)
pr = pr[pr$Organism=="Human",]

pr = data.frame(separate_rows(pr, "HGNC.ID", sep = "\\|"))
pr = data.frame(separate_rows(pr, "Entrez.Gene.ID", sep = "\\|"))
pr = data.frame(separate_rows(pr, "Ensembl.Gene.ID", sep = "\\|"))
pr = pr[!duplicated(pr),]

#obtaining the TSS information from Biomart in form of dataframe
t2g <- getBM(attributes=c("chromosome_name", "entrezgene_id", "hgnc_id", "hgnc_symbol", "transcript_is_canonical", "transcript_start", "transcript_end", "strand", "ensembl_gene_id","gene_biotype", "ensembl_transcript_id"),
      filters ="biotype",
      values  =c("protein_coding"),
      mart    = ensembl)

# keep only chroms  c(1:22,"X", "Y")
t2g = t2g[t2g$chromosome_name %in%  c(1:22,"X", "Y"),]
# Keep only canonical transcripts from Ensembl or longest transcript?
# t2g = t2g[!is.na(t2g$transcript_is_canonical),]

table(pr$HGNC.ID %in% t2g$hgnc_id)
table(pr$Entrez.Gene.ID %in% t2g$entrezgene_id)
table(pr$Ensembl.Gene.ID %in% t2g$ensembl_gene_id)

# ok, so we just make a decision then, if we go for most mapping, we will have to make the decision on what to do re. the multiple TSS and multiple genes per SomaID
# e.g. 
t2g[t2g$hgnc_symbol=="ZSCAN2",]

# write.table(my_ids.version, file="/group/diangelantonio/users/claudia/ensembl_SomaScan_V4.1_7K_Annotated_Content_20210616.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)






######################################
############# DECIDE WHAT TO DO BEFORE CONTINUING 
#######################################

# for matrixeqtl:
# geneid        chr     s1      s2
geneloc = pr2
names(geneloc) = c("geneid", "chr", "s1", "s2")
geneloc$chr = gsub("X", "22", geneloc$chr)
geneloc =geneloc[order(geneloc$chr, geneloc$s1),]
geneloc$chr = paste("chr", geneloc$chr, sep="")

write.table(geneloc, file="/group/diangelantonio/users/claudia/geneloc.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)


### try by entrez?
attributes$name[grep("entrez",attributes$name)]
t2g<-getBM(attributes=c('entrezgene_id','chromosome_name','start_position','end_position'), mart = ensembl)
my_ids <- data.frame(entrezgene_id=pr$Entrez.Gene.ID)
# my_ids$ensembl_gene_id <- gsub("\\..*","", my_ids$ensembl_gene_id_version)
my_ids.version <- merge(my_ids, t2g, by= 'entrezgene_id') # 7270 matching by entrezgene_id

### try by gne nmae
attributes$name[grep("hgnc",attributes$name)]
t2g<-getBM(attributes=c('hgnc_id','chromosome_name','start_position','end_position'), mart = ensembl)
my_ids <- data.frame(hgnc_id=pr$HGNC.ID)
# my_ids$ensembl_gene_id <- gsub("\\..*","", my_ids$ensembl_gene_id_version)
my_ids.version <- merge(my_ids, t2g, by= 'hgnc_id') # 2923 matching by hgnc_symbol

genes <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), mart = ensembl)


############################# using gencode v43
library(rtracklayer)
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz
gene_file = "/group/diangelantonio/users/claudia/gencode.v43.basic.annotation.gtf.gz"

gtf <- rtracklayer::import(gene_file)
gtf_df=as.data.frame(gtf)
gtf_df =gtf_df[gtf_df$type =="gene",]
gtf_df$ensembl_gene_id = gsub("\\..*","", gtf_df$gene_id)
table(pr$Ensembl.Gene.ID %in% gtf_df$ensembl_gene_id) # 6740

#geneInfo$gene_name = gtf_df[match(geneInfo$geneid, gtf_df$gene_id), "gene_name"]
