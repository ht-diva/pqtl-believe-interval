#!/home/mfilosi/opt/miniconda3/envs/ipyenv/bin/python
from pygwas.gwas import GWAS
from pygwas.geneticdata import GenDataList
from pygwas.conditionalanalysis import ConditionalAnalysis
import numpy as np
import os


def main(p, resdir, imppanel='WES13K', ncores=1,
         genetic_conf='genetic_data.json', dryrun=False,
         prefix="", outpath=""):
    # Read in the GWAS results
    mygwas = GWAS(p=p, resdir=resdir,
                  load_sumstat=True)

    # Get genetic data
    gd = GenDataList.from_json(genetic_conf)

    # Check hits in the tophits file
    hits = (mygwas.tophits['LOG10P'] > -np.log10(5e-8)).sum()
    if hits == 0:
        print("No hits found")
        exit()
    else:
        print(f"Running conditional analysis on {ncores} cores!")

    # Create and start the Conditional analysis class
    # - it will create the needed command to run separately based on
    #   chromosomes and phenotypes
    cc = ConditionalAnalysis(prefix=prefix)
    cc.create_cmd(gwas=mygwas, genetic_data=gd.gendata[imppanel],
                  outdir=outpath)

    # Run command for conditiaonl analysis (GCTA)
    cc.run(ncores=ncores, dryrun=dryrun)
    print(f"Run {len(cc.phenos)} phenotypes, merge files and cleanup logs...")
    cc.summarize_by_pheno()

    # Remove log files for each chromosome
    rem_files = cc.cleanup()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("phenotype")
    parser.add_argument("resdir")
    parser.add_argument("-c", "--cores", default=1, type=int)
    parser.add_argument("-g", "--genotype", default="WES13K")
    parser.add_argument("-d", "--dryrun", default=False, type=bool)
    parser.add_argument("--conf", default="genetic_data.json")
    parser.add_argument('-o', "--outpath", default="conditional_analysis")

    args = parser.parse_args()
    main(args.phenotype, args.resdir, imppanel=args.genotype,
         ncores=args.cores,
         genetic_conf=args.conf, dryrun=args.dryrun, outpath=args.outpath)
