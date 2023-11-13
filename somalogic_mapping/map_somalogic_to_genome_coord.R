# library(ensembldb)
library(dplyr)
library(tidyfr)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Hsapiens.v75)
library(glue)

#---- Setup ----#
proj_path <- "."
data_path <- glue("{proj_path}/data")

#---- EnsemblDB init ----
ens38 <- EnsDb.Hsapiens.v86 # GRCh38
ens37 <- EnsDb.Hsapiens.v75 # GRCh37

#---- Get somalogic metadata ----
ll <- read.table("somalogic_map_chris.csv", sep="\t", 
    header=TRUE) 

#---- Get uniprot ID to match ----#
mygeneid <- ll %>% dplyr::filter(UniProt != "") %>% dplyr::pull(EntrezGeneSymbol)
myuniprotid <- ll %>% dplyr::filter(UniProt != "") %>% dplyr::pull(UniProt)

#---- Match proteinID to Ensembl ----
myprot <- genes(ens37, filter=GeneNameFilter(mygeneid), columns = c(listProteinColumns(ens37), "symbol"))
mcols(myprot)$start <- start(myprot)
mcols(myprot)$end <- end(myprot)
mcols(myprot)$CHROM <- seqnames(myprot)

#---- Aggregate and find genomic coordinate for the protein/genes ----
genepos_df <- mcols(myprot) %>% as_tibble()

#---- Handle Multiple matches ----
# There are 2 genes that map to multiple position on the genome (different chromosomes)
# and have different ENSG ID. So they have to be treated separately in some way.
# I will assign select a specific ENSID for thos genes
# These are the genes with multiple matches on the ENSdb
genepos_df %>% 
  dplyr::filter(CHROM %in% as.character(c(1:22, "X", "Y"))) %>% 
  group_by(symbol) %>% 
  summarise(nid=n_distinct(gene_id), nchr=n_distinct(CHROM)) %>% 
  dplyr::filter(nid > 1, nchr>1)
multi_gene_id <- c(CKS1B="ENSG00000173207", LSP1="ENSG00000130592")

# Summarise by symbol, gene_id and CHROM and get max span of the gene
genepos_df <- genepos_df %>% dplyr::group_by(symbol, CHROM, gene_id) %>% 
  dplyr::summarise(start=min(start), end=max(end),
                   protein_id_n=n_distinct(protein_id), tx_id_n=n_distinct(tx_id)) %>% 
  ungroup() %>% 
  dplyr::filter(CHROM %in% as.character(c(1:22, "X", "Y")))
  
# Match multiple genes and index to remove from the DF
ix <- which(genepos_df$symbol %in% names(multi_gene_id) & !(genepos_df$gene_id %in% multi_gene_id))

# REMOVE unwanted hits
genepos_df <- genepos_df[-ix,]

# MERGE multiple hits
genepos_df <- genepos_df %>% group_by(symbol, CHROM) %>% summarise(start=min(start), end=max(end))

#---- Protein Position ----
prot_pos <- genepos_df %>% right_join(ll, by=c("symbol"="EntrezGeneSymbol"))

#---- Non matched proteins ----
# Get gene ids for each protein. If multiples gene id are present
# search for them separately and merge in a second time in ordert to maximize
# the overlap.
geneid_miss <- prot_pos %>% ungroup() %>% dplyr::filter(is.na(start), symbol != "") %>% 
  dplyr::select(label, EntrezGeneID) %>% 
  mutate(GIDL=str_split(EntrezGeneID, pattern="(\\|| )")) %>% 
  unnest(GIDL)

# Filter
mypromiss <- genes(ens37, filter=EntrezFilter(geneid_miss$GIDL), columns = c(listProteinColumns(ens37), "symbol"))
mcols(mypromiss)$start <- start(mypromiss)
mcols(mypromiss)$end <- end(mypromiss)
mcols(mypromiss)$CHROM <- seqnames(mypromiss)

#---- Aggregate and find genomic coordinate for the protein/genes ----
mypromiss <- mcols(mypromiss) %>% as_tibble()
mypromiss2 <- mypromiss %>% unnest(entrezid) %>% dplyr::group_by(symbol, CHROM, entrezid) %>% 
  dplyr::filter(CHROM %in% as.character(c(1:22, "X", "Y"))) %>% 
  dplyr::summarise(start=min(start), end=max(end),
                   protein_id_n=n_distinct(protein_id), tx_id_n=n_distinct(tx_id)) %>% 
  ungroup() %>% 
  mutate(entrezid=as.character(entrezid))

#---- Merge with original protein labels ----
prot_pos_miss <- geneid_miss %>% 
  left_join(mypromiss2, by=c("GIDL"="entrezid")) %>% 
  group_by(label, EntrezGeneID, CHROM, symbol) %>% 
  summarise(start=min(start), end=max(end)) %>%
  right_join(ll %>% 
               dplyr::filter(label %in% unique(geneid_miss$label)), 
             by=c("label", "EntrezGeneID"))

# Bind protein annotation from first merge attempt
prot_pos_max_intercept <- bind_rows(
  prot_pos %>% filter(!is.na(start)), 
  prot_pos_miss) 

# Save files
prot_pos %>% rename(start_grch37=start, end_grch37=end) %>% 
  write.table(file=file.path(proj_path, "somalogic_analysis", "protein_position.csv"),
              col.names = TRUE, row.names = TRUE, quote=FALSE, sep="\t")
prot_pos %>% rename(start_grch37=start, end_grch37=end) %>% 
  saveRDS(file=file.path(proj_path, "somalogic_analysis", "protein_position.rds"))

# Save file with the maximum intercept!
prot_pos_max_intercept %>% dplyr::rename(start_grch37=start, end_grch37=end) %>% 
  write.table(file=file.path(proj_path, "somalogic_analysis", "protein_position_max_intercept.csv"),
              col.names = TRUE, row.names = TRUE, quote=FALSE, sep="\t")
prot_pos_max_intercept %>% dplyr::rename(start_grch37=start, end_grch37=end) %>% 
  saveRDS(file=file.path(proj_path, "somalogic_analysis", "protein_position_max_intercept.rds"))

