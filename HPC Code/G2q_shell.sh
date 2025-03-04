#!/bin/bash

#PBS -N G2q
#PBS -M lucas.ostrowski@griffithuni.edu.au
#### select resources
#PBS -l walltime=999:00:00
#PBS -l select=1:ncpus=1:mem=24g
##### redirect error and output files
#PBS -e /export/home/s5257291/large_p_2024/G2q.err
#PBS -o /export/home/s5257291/large_p_2024/G2q.out

#### mail options
#### load matlab module (setup environment)
module load matlab/2021a
#### change to current working directory
cd /export/home/s5257291/large_p_2024
##### execute program
matlab -nodisplay -nodesktop -nosplash -r "G2_General(501, 50, 0, -1, 4*2*4*250^2/(pi^2*50), 101);" run.log

exit 0;