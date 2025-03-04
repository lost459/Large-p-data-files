#!/bin/bash

#PBS -N G1_lam_short
#PBS -M lucas.ostrowski@griffithuni.edu.au
#### select resources
#PBS -l walltime=999:00:00
#PBS -l select=1:ncpus=1:mem=20g
##### redirect error and output files
#PBS -e /export/home/s5257291/large_p_2024/G1_lam_short.err
#PBS -o /export/home/s5257291/large_p_2024/G1_lam_short.out

#### mail options
#### load matlab module (setup environment)
module load matlab/2021a
#### change to current working directory
cd /export/home/s5257291/large_p_2024
##### execute program
matlab -nodisplay -nodesktop -nosplash -r "G1_General(501, 50, 0.5, 0, 4*250^2/(pi^2*50), 120);" run.log

exit 0;