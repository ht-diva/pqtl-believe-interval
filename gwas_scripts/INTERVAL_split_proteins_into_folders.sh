###########################
### Claudia Giambartolomei
# 5 December 2023
###########################
# This nextflow pipeline was built by Edoardo Giacopuzzin the Genomics Lab
# https://github.com/HTGenomeAnalysisUnit/nf-pipeline-regenie
# set as false /exchange/healthds/pQTL/regenie/nf-pipeline-regenie/nextflow.config n_top_loci_plot

# What this script does:
# chunks phenotype file and creates folders for each chunk
# creates submission files adapted for each chunk
###########################
# INPUT
PhenoFile=/exchange/healthds/pQTL/INTERVAL/INTERVAL_NonImp_residuals_new_ReorderCols_FID_IID_OtherPheno.txt
main_conf=/exchange/healthds/pQTL/regenie/claudia_test/interval_dec_tests/split_proteins/single_project_targets_1_2.conf
main_sbatch=/exchange/healthds/pQTL/regenie/claudia_test/interval_dec_tests/split_proteins/regenie_test_targets.sbatch
main_dir=/exchange/healthds/pQTL/regenie/claudia_test/interval_dec_tests/split_proteins/
###########################
submit=false
###########################
cp $PhenoFile .
myfile=INTERVAL_NonImp_residuals_new_ReorderCols_FID_IID_OtherPheno.txt
chunksize=100
ncol=7148
# chunks of 100 are the most ideal
for ((col=3; col <= ncol ; col = col+chunksize))
do
        mkdir chunk_$col
        cut -f1,2,$col-$((col+chunksize-1)) < $myfile > chunk_$col/INTERVAL_NonImp_residuals_new_ReorderCols_FID_IID_OtherPheno.txt
done

chunks=$(ls -d chunk_*)
for i in $chunks; do
        echo $i;

        cd $main_dir
        cd chunk_$i
        chunk_conf=chunk_${i}.conf
        chunk_sbatch=chunk_${i}.sbatch

        cp $main_conf $chunk_conf
        cp $main_sbatch $chunk_sbatch
        pheno=$PWD/INTERVAL_NonImp_residuals_new_ReorderCols_FID_IID_OtherPheno.txt
        # print phenotypes comma separated

        sed 's/test_targets_1_2/chunk_'$i'/' $chunk_conf > temp && mv temp $chunk_conf
        sed 's#/exchange/healthds/pQTL/INTERVAL/INTERVAL_NonImp_residuals_new_ReorderCols_FID_IID_OtherPheno.txt#'$PWD'/INTERVAL_NonImp_residuals_new_ReorderCols_FID_IID_OtherPheno.txt#g' $chunk_conf > temp && mv temp $chunk_conf
        list_proteins=$(awk 'NR==1' $pheno | cut -f 3- | sed -e 's/\s\+/,/g')
        sed 's/seq.10000.28_res,seq.10001.7_res/'$list_proteins'/' $chunk_conf > temp && mv temp $chunk_conf

        sed 's/single_project_targets_1_2.conf/'$chunk_conf'/' $chunk_sbatch > temp && mv temp $chunk_sbatch
        sed 's/test_/test_chunk_'$i'/' $chunk_sbatch > temp && mv temp $chunk_sbatch
        
        if [ "$submit" = true ] ; then
                sbatch $chunk_sbatch
        fi

done
