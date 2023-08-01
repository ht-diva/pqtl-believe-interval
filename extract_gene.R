if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

# Aptamer level annotations were created by mapping proteins to genomic coordinates using GENCODE (GRCh38), version 32 (Ensembl 98) [76].
# GENCODE: http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&c=chrX&g=knownGene

pr_file = "/center/healthds/pQTL/Reference_datasets_for_QC_proteomics/Somalogic/SomaScan_V4.1_7K_Annotated_Content_20210616.csv"
pr = read.csv(pr_file, header=T,stringsAsFactors=F)

biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl))
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position'), mart = ensembl)
my_ids <- data.frame(ensembl_gene_id=pr$Ensembl.Gene.ID)
# my_ids$ensembl_gene_id <- gsub("\\..*","", my_ids$ensembl_gene_id_version)
my_ids.version <- merge(my_ids, t2g, by= 'ensembl_gene_id') # 6746 matching by ensembl_gene_id

pr2 = merge(pr, my_ids.version, by.x = "Ensembl.Gene.ID", by.y="ensembl_gene_id")
pr2 = pr2[,c("X...SeqId", "chromosome_name", "start_position", "end_position")]
pr2 = pr2[!duplicated(pr2),] # 6746

# write.table(my_ids.version, file="/group/diangelantonio/users/claudia/ensembl_SomaScan_V4.1_7K_Annotated_Content_20210616.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

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
