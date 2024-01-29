import gwaslab as gl
import argparse
import json
from multiprocessing import Pool
import numpy as np



def process_file(f, config):
    f = str(f)
    study = ".".join(f.split('.')[:-3][1:])
    
    print(f"Processing {f}")
    
    mysumstats = gl.Sumstats(f,
                             study=study,
                             snpid="ID",
                             chrom="CHROM",
                             pos="GENPOS",
                             ea="ALLELE0",
                             nea="ALLELE1",
                             eaf="A1FREQ",
                             beta="BETA",
                             se="SE",
                             mlog10p="LOG10P",
                             verbose=False)
        
    save_name = config['save_path'] + study + '.png'
    dist = config['cut_dist']
    cut = round(-np.log10(config['sig_level'])) + dist

    mysumstats.plot_mqq(skip=config['skip'],
                        cut=cut,
                        scaled=config['scaled'], 
                        build=config['build'], 
                        mode="mqq", 
                        stratified=config['stratified'], 
                        maf_bins=config['maf_bins'], 
                        figargs={"figsize":(15,5),"dpi":300},
                        marker_size=(1,7),
                        save=save_name,
                        sig_level=config['sig_level'],
                        verbose=False)

    
def main(config_path):
    
    with open(config_path) as config_file:
        config = json.load(config_file)

    with open(config['gwas_files_list_dir']) as file:
        file_list = [line.strip() for line in file]

    pool = Pool(processes=config['n_processes'])  # Adjust the number of processes as needed
    for f in file_list:
        pool.apply_async(process_file, args=(f, config))
        
    pool.close()
    pool.join()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process GWAS data.')
    parser.add_argument('config', type=str, help='Path to JSON config file')
    args = parser.parse_args()

    main(args.config)
