#!/bin/bash
##############################################################################
#BATCH --job-name=4_hayridurmaz-rowmv
#SBATCH --account=users
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --partition=cpu
#SBATCH --output=4_hayridurmaz-rowmv.out
##############################################################################
hostname
# Module load
module load openmpi/4.0.3

cd /home/users/hayridurmaz/rowmv
make
# Execute the program
srun -n 4 /home/users/hayridurmaz/rowmv/rowmv
