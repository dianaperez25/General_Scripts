#!/bin/bash 
#SBATCH -A b1081
#SBATCH -n 40
#SBATCH --mem=0
#SBATCH -t 50:00:00 
#SBATCH -p b1081
#SBATCH --job-name="LatPermutations"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dianaperezrivera2024@u.northwestern.edu
#SBATCH -o "%x.o%j"

## set your working directory 
cd $PATH:/projects/p31161/lateralizationVariants/

## job commands; <matlabscript> is your MATLAB .m file, specified without the .m extension
module load matlab/r2016a
matlab -nosplash -nodesktop -singleCompThread -r permute_diffMaps_flippedVars
