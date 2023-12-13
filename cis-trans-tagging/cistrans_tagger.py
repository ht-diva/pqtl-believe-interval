import argparse
import pandas as pd
import os
import glob
import json
import mpmath  # Make sure this module is installed
import logging
import warnings



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

    

def cis_trans_tagger(config_file):
    # Load configuration
    params = parse_config(config_file)
    
    # Set up logging
    log_file_path = os.path.join(params['save_path'], 'cis_trans_tagging_log.txt')
    setup_logging(log_file_path)

    # Load mapping file 
    map_df = pd.read_csv(params['map_file_path'], index_col='Unnamed: 0')
    seqid_list = map_df[params['id_column_name']].unique()

    # column types needed for safe loading of specific columns
    column_types = {params['snp_file_chrom_col']: str,
                    params['pval_column']: str}

    # Initialize an empty DataFrame for appending results
    all_snps_combined_df = pd.DataFrame()

    # Process each target protein
    for i, seqid in enumerate(seqid_list, start=1):
        
        print(f"Processing SeqId {seqid} ({i} of {len(seqid_list)})", end='\r')

        target_protein = map_df[map_df[params['id_column_name']] == seqid]
        if not target_protein.empty:
            target_chrom = target_protein.iloc[0][params['map_file_chrom_col']]
            start_position = target_protein.iloc[0][params['map_file_start_col']]
            end_position = target_protein.iloc[0][params['map_file_end_col']] + params['region_buffer_value']

            # Find pQTL results files
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
                    snp_df.loc

                    # Calculate distance from bounds and identify proximal signals
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

                    snp_df[params['id_column_name']] = seqid

                    # Append to the combined DataFrame
                    all_snps_combined_df = pd.concat([all_snps_combined_df, snp_df], axis=0)
                    break
            
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
