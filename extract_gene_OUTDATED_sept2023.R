# super useful function 'separate_rows' from tidyr to separate rows by delimiter
library(tidyr)
#pr = data.frame(separate_rows(pr, "HGNC.ID", sep = "\\|"))
#pr = data.frame(separate_rows(pr, "Entrez.Gene.ID", sep = "\\|"))
#pr = data.frame(separate_rows(pr, "Entrez.Gene.Name", sep = "\\|"))
#pr = data.frame(separate_rows(pr, "Ensembl.Gene.ID", sep = "\\|"))
#pr = data.frame(separate_rows(pr, "UniProt.ID", sep = "\\|"))

#############################
# Mapped using Entrez gene ID
# 7289 SeqIDs and 7377 SeqID_EntrezGeneName combo
pr_file = "/center/healthds/pQTL/Reference_datasets_for_QC_proteomics/Somalogic/SomaScan_V4.1_7K_Annotated_Content_20210616.csv"
pr = read.csv(pr_file, header=T,stringsAsFactors=F)
table(pr$Organism)
pr = pr[pr$Organism=="Human",]
# 20367-6, remove b/c not human protein (Type)
# super useful function 'separate_rows' from tidyr to separate rows by delimiter
pr = pr[,c("X...SeqId", "SomaId", "Target.Name", "Target.Full.Name", "UniProt.ID", "Entrez.Gene.Name", "Entrez.Gene.ID", "Ensembl.Gene.ID", "HGNC.ID")]

############################
# We used biomaRt (Ensembl Genes 110, dataset Human genes version GRCh38.p14)
# We first matched the SomaScan HGNC.ID (for 7289 SeqIDs and 6378 HGNC.ID) 
# For the 6 SeqIDs we could not find a match using HGNC.ID, we proceeded subsequently looking for match using Entrez.Gene.ID.
# For 2, we could not find any matches for TSS. 

# a match between hgnc_symbol and Entrez.Gene.Name. 
# for the 15 SeqIDs we could still not find a match using HGNC.ID, we proceeded subsequently looking for a match using:
# 2.

	library(biomaRt)
	ensembl <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
	# obtaining the TSS information from Biomart in form of dataframe
	t2g_biomart = getBM(attributes=c("chromosome_name", "entrezgene_id", "hgnc_id", "hgnc_symbol", "transcript_is_canonical", "transcription_start_site", "strand", "ensembl_gene_id","gene_biotype"),
      	mart    = ensembl)
	t2g = t2g_biomart[!is.na(t2g_biomart$transcript_is_canonical),]
	t2g = t2g[t2g$chromosome_name %in% c(1:22,"X", "Y"),]
	t2g = t2g[t2g$gene_biotype == "protein_coding",]

	###
	pr_hgnc = data.frame(separate_rows(pr, "HGNC.ID", sep = "\\|"))
	pr_hgnc = pr_hgnc[pr_hgnc$HGNC.ID !="",]

	matched_hgnc = merge(pr_hgnc, t2g, by.x = "HGNC.ID", by.y = "hgnc_id")
	matched_hgnc$match = "HGNC.ID"
                table(matched_hgnc$chromosome_name)
                table(matched_hgnc$gene_biotype)

	#not_matched = pr_hgnc[!(pr_hgnc$HGNC.ID %in% t2g$hgnc_id),]
	not_matched = pr[!(pr$X...SeqId %in% matched_hgnc$X...SeqId),]

  	print(length(unique(not_matched$X...SeqId)))

        ###
        pr_entrez = data.frame(separate_rows(not_matched, "Entrez.Gene.ID", sep = "\\|"))
	pr_entrez = pr_entrez[pr_entrez$Entrez.Gene.ID !="",]
        matched_entrez = merge(pr_entrez, t2g, by.x = "Entrez.Gene.ID", by.y = "entrezgene_id")
	matched_entrez$match = "Entrez.Gene.ID"
                table(matched_entrez$chromosome_name)
                table(matched_entrez$gene_biotype)
	not_matched = pr_entrez[!(pr_entrez$X...SeqId %in% matched_entrez$X...SeqId),]
        #not_matched = pr_entrez[!(pr_entrez$Entrez.Gene.ID %in% t2g$entrezgene_id),]
                print(length(unique(not_matched$X...SeqId)))

	### This does not add anything 
	go = FALSE
	if (go) {
	pr_entrezName = data.frame(separate_rows(not_matched, "Entrez.Gene.Name", sep = "\\|"))
	pr_entrezName = pr_entrezName[pr_entrezName$Entrez.Gene.Name != "",]
	matched_entrezName = merge(pr_entrezName, t2g, by.x = "Entrez.Gene.Name", by.y = "hgnc_symbol")
	matched_entrezName$match= "Entrez.Gene.Name"
  		table(matched_entrezName$chromosome_name)
  		table(matched_entrezName$gene_biotype)
	not_matched = pr_entrezName[!(pr_entrezName$X...SeqId %in% matched_entrezName$X...SeqId),]
	#not_matched = matched_entrezName[!(matched_entrezName$Entrez.Gene.Name %in% t2g$hgnc_symbol),]
  		print(length(unique(not_matched$X...SeqId)))
	}

	###
        go = FALSE
        if (go) {
	pr_ensembl = data.frame(separate_rows(not_matched, "Ensembl.Gene.ID", sep = "\\|"))
	pr_ensembl = pr_ensembl[pr_ensembl$Ensembl.Gene.ID != "",]
	matched_ensembl = merge(pr_ensembl, t2g, by.x = "Ensembl.Gene.ID", by.y = "ensembl_gene_id")
	matched_ensembl$match = "Ensembl.Gene.ID"
  		table(matched_ensembl$chromosome_name)
  		table(matched_ensembl$gene_biotype)
	not_matched = pr_ensembl[!(pr_ensembl$X...SeqId %in% matched_ensembl$X...SeqId),]
	#not_matched = pr_ensembl[!(pr_ensembl$Ensembl.Gene.ID %in% t2g$ensembl_gene_id),]
  		print(length(unique(not_matched$X...SeqId)))
	}

#############
# MERGE THE MATCHED DATA

cols = c(names(pr), "match", "chromosome_name", "strand", "transcription_start_site","gene_biotype")
matched_hgnc = matched_hgnc[,cols]
matched_entrez = matched_entrez[,cols]
matched = rbind.data.frame(matched_hgnc,matched_entrez)

not_matched$match = "lookup"
not_matched$chromosome_name = NA
not_matched$transcription_start_site = NA
not_matched$gene_biotype = NA

#matched = matched[matched$chromosome_name %in% c(1:22,"X", "Y"),]
length(unique(matched$X...SeqId))
# 7254 SeqIDs and 35 SeqIDs no match
# Use SeqId_Entrez.Gene.Name as unique combination
matched$SeqId_Entrez.Gene.Name = paste(matched$X...SeqId, matched$Entrez.Gene.Name, sep="_")
matched = matched[!duplicated(matched),]

not_matched = pr[!(pr$X...SeqId %in% matched$X...SeqId),]
not_matched$SeqId_Entrez.Gene.Name = paste(not_matched$X...SeqId, not_matched$Entrez.Gene.Name, sep="_")
# import the info from other file to see if i can recover some?
## gene_file =  "/group/diangelantonio/users/claudia/biomaRt_comparison_SomaScan_V4.1_7K_Annotated_Content_20210616.txt"
## gene=read.delim(gene_file, header=T, stringsAsFactors=F, sep="\t")
#f = merge(not_matched, gene, by = "SeqId_Entrez.Gene.Name")

#not_matched$chromosome_name[not_matched$Entrez.Gene.Name=="TXNRD3NB"] = 
#?
#not_matched$chromosome_name[not_matched$Entrez.Gene.Name=="BAGE3"] = chr12:111471309
#TXNRD3NB"chr3:126,571,779


#geneid  chr     s1      s2
#/group/diangelantonio/users/claudia/geneloc.txt


write.table(matched, file="/exchange/healthds/pQTL/pQTL_workplace/matched_7256_biomaRt_SomaScan_V4.1_7K_Annotated_Content_20210616.txt", sep="\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

