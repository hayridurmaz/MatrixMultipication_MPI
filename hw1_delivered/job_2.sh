#!/bin/bash
##############################################################################
#BATCH --job-name=2_hayridurmaz-rowmv
#SBATCH --account=users
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cpu
#SBATCH --output=2_hayridurmaz-rowmv.out
##############################################################################
hostname
# Module load
module load openmpi/4.0.3

cd /home/users/hayridurmaz/rowmv
make
# Execute the program
srun -n 2 /home/users/hayridurmaz/rowmv/rowmv
