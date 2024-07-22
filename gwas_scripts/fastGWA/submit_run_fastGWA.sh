main=run_fastGWA.sh

  
wk_dir=temp_scripts
mkdir -p $wk_dir
cd $wk_dir

#for colNum in {3..7246}; do
#for colNum in {3..9219}; do # a lot more!
#for colNum in {3..100}; do

for colNum in {3..368}; do
# 3:368 tests! (366 tests!)


#for colNum in {111..1990}; do
#for colNum in {1991..3000}; do
#for colNum in {3000..4900}; do
#for colNum in {4900..6000}; do
#for colNum in {6001..7246}; do
        script=${colNum}.cmd
        cp $main $script
        #sed 's/chr=1/chr='$chr'/' $script > temp && mv temp $script
        sed 's/colNum=3/colNum='$colNum'/' $script > temp && mv temp $script
        #sed 's/PBS -N coloc/PBS -N coloc_'$n_comb'/' $script > temp && mv temp $script

        sbatch $script
        #qsub -t 1-$njobs:1 $script
done

