import argparse
import pandas as pd
import os
import glob
import json
import mpmath  
import logging
import warnings
import pybedtools



# Function to transform strings to numbers using mpmath
def to_mpf(x):
    return mpmath.mpf(x)

# Function to read and parse the config file
def parse_config(config_file):
    with open(config_file, 'r') as file:
        return json.load(file)
    
# Clear screen function
def clear_screen():
    os.system('cls' if os.name == 'nt' else 'clear')
    
    
# Function to configure logging
def setup_logging(log_file_path):
    logging.basicConfig(filename=log_file_path, level=logging.INFO,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # Capture warnings with logging
    warnings.showwarning = lambda message, category, filename, lineno, file=None, line=None: \
        logging.warning('%s:%s: %s:%s', filename, lineno, category.__name__, message)

    
    
def get_unique_trans_snps(all_snps_df, params):
    # Extract unique trans SNPs based on chromosome and position
    unique_trans_snps = all_snps_df[all_snps_df['cis_trans'] == 'trans'].drop_duplicates(subset=[params['snp_file_chrom_col'], params['snp_file_pos_col']])
    return unique_trans_snps

       
def batch_annotate_trans_snps(trans_snps_df, genes_bed, params):
    # Create a BED file from trans SNPs DataFrame
    trans_snps_bed = pybedtools.BedTool.from_dataframe(trans_snps_df[[params['snp_file_chrom_col'], params['snp_file_pos_col'], params['snp_file_pos_col']]])

    # Annotate with nearest gene in batch
    closest_genes = trans_snps_bed.closest(genes_bed, d=True).to_dataframe()
    
    return closest_genes[param['gene_file_colname']]

# helper function to create BED file from DataFrame for pybedtools
def df_to_bed(df, chrom_col, start_col, end_col, name_col=None):
    bed_df = df[[chrom_col, start_col, end_col]]
    if name_col:
        bed_df = bed_df.assign(name=df[name_col])
    return pybedtools.BedTool.from_dataframe(bed_df)



def cis_trans_tagger(config_file):
    ## Load configuration
    params = parse_config(config_file)

    # Set up logging
    log_file_path = os.path.join(params['save_path'], 'cis_trans_tagging_log.txt')
    setup_logging(log_file_path)


    # Load mapping file with pandas
    map_df = pd.read_csv(params['target_gene_map_file_dir'], index_col='Unnamed: 0')
    seqid_list = map_df[params['prot_id_column_name']].unique()

    # column types needed for safe loading of results
    column_types = {params['snp_file_chrom_col']: str,
                    params['pval_column']: str}

    # Initialize an empty DataFrame for appending results
    all_snps_combined_df = pd.DataFrame()

    # Process each SeqId
    for i, seqid in enumerate(seqid_list, start=1):

        print(f"Processing SeqId {seqid} ({i} of {len(seqid_list)})", end='\r')

        target_protein = map_df[map_df[params['prot_id_column_name']] == seqid]
        if not target_protein.empty:
            target_chrom = target_protein.iloc[0][params['map_file_chrom_col']]
            start_position = target_protein.iloc[0][params['map_file_start_col']]
            end_position = target_protein.iloc[0][params['map_file_end_col']] + params['region_buffer_value']

            # Find SNP files
            snp_files = glob.glob(os.path.join(params['gwas_files_dir'], f'*{seqid}*'))

            # Check for multiple files
            if len(snp_files) > 2:
                clear_screen()
                print(f"Too many files found for SeqId {seqid}. Check logs for details. Processing first file only.")
                logging.warning(f"Too many files found for SeqId {seqid}. Processing the first encountered file only. Files found: {', '.join(snp_files)}")

            for file in snp_files:
                if '.log' not in file:
                    snp_df = pd.read_csv(file, sep="\t", dtype=column_types)

                    # transform strings to numbers
                    snp_df[params['pval_column']] = snp_df[params['pval_column']].apply(to_mpf)

                    # Apply P-value threshold filter
                    snp_df = snp_df[snp_df[params['pval_column']] < params['pvalue_threshold']]

                    # Determine cis/trans
                    on_same_chrom = snp_df[params['snp_file_chrom_col']] == target_chrom
                    in_range = snp_df[params['snp_file_pos_col']].between(start_position, end_position) 
                    snp_df['cis_trans'] = 'trans'
                    snp_df.loc[on_same_chrom & in_range, 'cis_trans'] = 'cis'

                    snp_df['nearest_gene'] = pd.NA
                    cis_snps = snp_df[snp_df['cis_trans'] == 'cis']
                    for index, row in cis_snps.iterrows():
                        gene_name = target_protein.iloc[0][params['map_file_gene_name']]
                        snp_df.at[index, 'nearest_gene'] = gene_name

                    # Calculate distance from bounds
                    snp_df['dist_lower'] = pd.NA
                    snp_df['dist_upper'] = pd.NA
                    snp_df.loc[on_same_chrom, 'dist_lower'] = abs(snp_df[params['snp_file_pos_col']] - start_position)
                    snp_df.loc[on_same_chrom, 'dist_upper'] = abs(snp_df[params['snp_file_pos_col']] - end_position)
                    proximal = (snp_df['dist_lower'] < params['proximality_threshold']) | (snp_df['dist_upper'] < params['proximality_threshold'])
                    snp_df['proximal'] = 'N'
                    snp_df.loc[proximal, 'proximal'] = 'Y'

                    # Drop unnecessary columns
                    columns_to_keep = set(snp_df.columns) - set(params['columns_to_drop'])
                    snp_df = snp_df[list(columns_to_keep)]

                    snp_df[params['prot_id_column_name']] = seqid

                    # Append to the combined DataFrame
                    all_snps_combined_df = pd.concat([all_snps_combined_df, snp_df], axis=0)
                    break


    print("\nMoving to annotating trans.")                

    # Load the data into a DataFrame
    gene_df = pd.read_csv(params['gene_pos_file_dir'], sep='\s+', header=0)

    # Ensuring chromosome is a string for sorting
    gene_df['CHR'] = gene_df['CHR'].astype(str)

    # Step 1: Filter Unique Trans SNPs
    unique_trans_snps = all_snps_combined_df[all_snps_combined_df['cis_trans'] == 'trans'].drop_duplicates(subset=[params['snp_file_chrom_col'], params['snp_file_pos_col']])
    unique_trans_snps[params['snp_file_chrom_col']] = unique_trans_snps[params['snp_file_chrom_col']].astype(str)

    # Sorting the DataFrames
    unique_trans_snps = unique_trans_snps.sort_values(by=[params['snp_file_chrom_col'], params['snp_file_pos_col']])
    gene_df = gene_df.sort_values(by=['CHR', 'START'])


    # Create BED file for SNPs
    snps_bed = df_to_bed(unique_trans_snps, params['snp_file_chrom_col'], params['snp_file_pos_col'], params['snp_file_pos_col'])

    # Create BED file for genes, including gene names
    genes_bed = df_to_bed(gene_df, 'CHR', 'START', 'END', 'GENE_NAME')


    # Use pybedtools to find the closest genes
    closest_genes = snps_bed.closest(genes_bed, d=True, t='first').to_dataframe(names=['chr', 'start', 'end', 'gene_chr', 'gene_start', 
                                                                                       'gene_end', 'gene_name', 'score', 'strand', 'distance'])

    # Creating a dictionary for mapping with keys as tuples (chromosome, position)
    snp_to_gene_map = dict(zip(zip(closest_genes['chr'].astype(str), closest_genes['start']), closest_genes['gene_name']))

    # Update the mapping logic to correctly reference separate chromosome and position columns
    for index, row in all_snps_combined_df.iterrows():
        if row['cis_trans'] == 'trans':
            chrom = str(row[params['snp_file_chrom_col']])
            pos = row[params['snp_file_pos_col']]
            if (chrom, pos) in snp_to_gene_map:
                all_snps_combined_df.at[index, 'nearest_gene'] = snp_to_gene_map[(chrom, pos)]



    initial_cols = [params['prot_id_column_name'],  params['snp_file_chrom_col'], params['snp_file_pos_col'], params['pval_column'], 'cis_trans', 'proximal', 'dist_lower', 'dist_upper', 'nearest_gene']
    other_columns = [col for col in all_snps_combined_df.columns if col not in initial_cols]
    # Combine the lists to get the new column order
    new_column_order = initial_cols + other_columns

    all_snps_combined_df = all_snps_combined_df[new_column_order]
    
    print("Processing completed. Saving results table.")
    
    # Save the DataFrame as a pickle file
    pickle_file_path = os.path.join(params['save_path'], params['output_filename'])
    all_snps_combined_df.to_pickle(pickle_file_path)
    
    print("All done.")


    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Tag cis or trans pQTL signals from clumped results files.")
    parser.add_argument("config_file", help="Path to the configuration file")
    args = parser.parse_args()

    cis_trans_tagger(args.config_file)
