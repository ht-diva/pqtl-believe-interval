#########################
# Claudia Giambartolomei 
# 13 June 2024
#########################
# Extract genotypes/dosages and phenotypes and make boxplots by geno group
#########################
# INPUT
#snp = "1:196655713:C:T"
#seqid = "seq.12756.3"
#path_genFile="/scratch/dariush.ghasemi/projects/genomics_QC_pipeline/results/updated_snpid/chr1"
# "/scratch/giulia.pontali/genomics_QC_pipeline/results/pgen/impute_recoded_selected_sample_filter_hq_var_3"
#path_phenFile = "/exchange/healthds/pQTL/results/INTERVAL/INTERVAL_NonImp_residuals_final.txt"

##########################
plot_snpid_seqid = function(snp = "1:196655713:C:T", seqid = "seq.12756.3", path_genFile="/scratch/dariush.ghasemi/projects/genomics_QC_pipeline/results/updated_snpid/chr1", path_phenFile = "/exchange/healthds/pQTL/results/INTERVAL/INTERVAL_NonImp_residuals_final.txt", path_plink="/home/c.giambartolomei/softwares/plink2.sh", out_dir="/scratch/c.giambartolomei/TEST_META/extraction/geno_plots/") {

		# load libraries
		library(ggplot2)
		library(plink2R)
		library(data.table)
		library(dplyr)

		if (!file.exists(out_dir)) dir.create(out_dir)

		prefix = paste(snp, seqid, sep="_")
                temp_prefix = paste(prefix, "_temp", sep="")
                temp_pfile = paste(out_dir, temp_prefix, sep="")
		temp_datafile = paste(out_dir, temp_prefix, ".txt", sep="")
		plot_file = paste(out_dir, prefix, ".png", sep="")

                if (!file.exists(paste(path_genFile, ".bim", sep=""))) stop(paste("File does not exist: ", path_genFile, ".bim", sep=""))

               
                cmd1 = paste(path_plink, " -bfile ", path_genFile, " --snp ", snp, " --make-bed --out ", temp_pfile, sep="")
		cmd1_error  <- tryCatch(system(cmd1), error=function(e) NULL )

		if (cmd1_error!=0) stop("The genotypes could not be extracted")
	
		# import geno	
		library(plink2R)
		genos = read_plink(temp_pfile)
		geno_data = as.data.frame(genos$bed)
		snpid = colnames(geno_data)
		rownames(geno_data) = genos$fam$V1
		sampleid_geno = rownames(geno_data)
		# import pheno
		pheno_names = scan(path_phenFile, what="character", nlines=1)
		ind = which(pheno_names == seqid)
		message("Reading phenotype data...")
		pheno_data = fread(path_phenFile, header=T, stringsAsFactors=F, data.table=F, select = c(1,2,ind))
		# check you imported the correct column and rename
		if (colnames(pheno_data)[3] != seqid ) stop("Import of pheno file did not work: seqid not in pheno file?")
		colnames(pheno_data)[3] = "pheno"
		pheno_data$IID = as.character(pheno_data$IID)
		sampleid_pheno = pheno_data$IID
		
		if (!all(sampleid_geno == sampleid_pheno)) stop("sample names don't match")

		# Merge phenotype data with genotype data
		pheno_data$geno_groups = geno_data[match(pheno_data$IID, rownames(geno_data)),]
		# recode geno groups into increasing levels of EFFECT allele
		EA = unlist(strsplit(snpid, ":"))[3]
		NEA = unlist(strsplit(snpid, ":"))[4]
		hom_minor = paste(EA, EA, sep="")
		het = paste(EA, NEA, sep="")
		hom_major = paste(NEA, NEA, sep="")

		pheno_data = pheno_data %>% mutate(geno_groups2=recode(geno_groups, '0'=hom_minor, '1'=het, '2'=hom_major)) %>% data.frame()
		pheno_data$geno_groups2[is.na(pheno_data$geno_groups2)] <- "NA"
		pheno_data$geno_groups2 <- factor(pheno_data$geno_groups2, levels = c(hom_minor, het, hom_major, "NA"))

		# save?
		message("File saved in: ", temp_datafile)
		write.table(x=pheno_data, file=temp_datafile, row.names = FALSE, col.names=TRUE, quote = FALSE, sep = '\t')

		# Calculate sample size for each genotype group
		sample_size <-	pheno_data %>%
  				group_by(geno_groups2) %>%
  				summarise(n = n())
		# Function to add sample size annotation
		add_sample_size <- function(p, data) {
  				sample_size <- data %>%
    				group_by(geno_groups2) %>%
    				summarise(n = n())

		   		p + geom_text(data = sample_size, aes(x = geno_groups2, y = max(data$pheno) * 1.05, label = paste("n =", n)), 
                		vjust = -0.5, size = 4)
		}

		# Create the plot
		p <- ggplot(pheno_data, aes(x = geno_groups2, y = pheno, color = geno_groups2)) +
  			geom_boxplot(outlier.shape = NA, fill = "transparent") +
  			geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  			stat_summary(fun = mean, geom = "point", shape = 4, size = 3, color = "red") +
  			labs(title = "Protein Levels by Genotype",
       			x = "Genotype",
       			y = "Protein Level") +
  			theme_minimal(base_size = 15) +
  			theme(legend.position = "none")

		# Add sample size annotation
		p <- add_sample_size(p, pheno_data)


		# Display the plot
		png(plot_file)
		print(p)
		dev.off()
		message("File in: ", plot_file)

		clean = T
		if (clean) {
		do.call(file.remove, list(list.files(out_dir, full.names = TRUE, pattern=temp_prefix)))
		}
	}
