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
    mygwas = GWAS(p=p, resdir=resdir, load_sumstat=True)#, imppanel=imppanel)
    # tpfile = os.path.join(mygwas.resdir, 'tophits', f'{mygwas.p}.regenie.filtered.gz')
    # mygwas.tophits = tpfile

    # Get genetic data
    gd = GenDataList.from_json(genetic_conf)

    # Check hits in the tophits file
    nhits = (mygwas.tophits['LOG10P'] > -np.log10(5e-8)).sum()
    if nhits == 0:
        print("Ziobilli")
        exit()
    else:
        print("Che cazzo.....")

    cc = ConditionalAnalysis(prefix=prefix)
    cc.create_cmd(gwas=mygwas, genetic_data=gd.gendata[imppanel],
                  outdir=outpath)

    print(ncores)

    # cc.run(ncores=ncores, dryrun=dryrun)
    print(f"Summarize for the pheno: {p}")
    # cc.summarize_by_pheno()
    rem_files = cc.cleanup()

    print(rem_files)
    print(rem_files['files'][rem_files['removed']].apply(os.path.basename))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("phenotype")
    parser.add_argument("resdir")
    parser.add_argument("-c", "--cores", default=1, type=int)
    parser.add_argument("-g", "--genotype", default="WES13K")
    parser.add_argument("-d", "--dryrun", default=False, type=bool)
    parser.add_argument("--conf", default="genetic_data.json")
    parser.add_argument('-o',  "--outpath", default="conditional_analysis")

    args = parser.parse_args()
    main(args.phenotype, args.resdir, imppanel=args.genotype, ncores=args.cores,
         genetic_conf=args.conf, dryrun=args.dryrun, outpath=args.outpath)
