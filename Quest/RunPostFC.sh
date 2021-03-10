#!/bin/bash 
#SBATCH -A b1081
#SBATCH -n 40
#SBATCH --mem=0
#SBATCH -t 20:00:00 
#SBATCH -p b1081
#SBATCH --job-name="TaskMatched"
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dianaperezrivera2024@u.northwestern.edu
#SBATCH -o "%x.o%j"

## set your working directory 
cd $PATH:/projects/p31161/Scripts/

## job commands; <matlabscript> is your MATLAB .m file, specified without the .m extension
module load matlab/r2016a
matlab -nosplash -nodesktop -singleCompThread -r run_postFC
