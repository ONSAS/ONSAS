#!/bin/bash
#SBATCH --job-name=staticVonMisses
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4096
#SBATCH --time=00:30:00
#SBATCH --tmp=9G
#SBATCH --partition=normal
#SBATCH --qos=normal

#PATH TO MATLAB bin: /clusteruy/apps/matlab/R2018b/bin/matlab  
#ALIAS FOR MATLAB bin: alias matlab = "/clusteruy/apps/matlab/R2018b/bin/matlab"
/clusteruy/apps/matlab/R2018b/bin/matlab -nodisplay -nosplash -nodesktop -r "run('./onsasExample_staticVonMisesTruss.m');exit;"