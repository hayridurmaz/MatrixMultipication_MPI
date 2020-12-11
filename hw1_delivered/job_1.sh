#!/bin/bash
##############################################################################
#BATCH --job-name=1_hayridurmaz-rowmv
#SBATCH --account=users
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cpu
#SBATCH --output=1_hayridurmaz-rowmv.out
##############################################################################
hostname
# Module load
module load openmpi/4.0.3

cd /home/users/hayridurmaz/rowmv
make
# Execute the program
srun -n 1 /home/users/hayridurmaz/rowmv/rowmv
