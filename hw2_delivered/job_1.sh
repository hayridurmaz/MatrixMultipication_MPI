#!/bin/bash
###################################################################################################
#SBATCH --job-name=1_hayridurmaz-rowmv
#SBATCH --account=users
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=cpu
#SBATCH --output=1_hayridurmaz-rowmv.out
###################################################################################################

module load gcc/9.2.0
module load openmpi/4.0.3
echo $HOSTNAME
echo $SLURM_CPUS_PER_TASK
make
#cat /proc/cpuinfo
for N in 1000 2000 4000 8000 16000 32000 64000; do
        for nthreads in 32; do
                export OMP_NUM_THREADS=$nthreads
                srun -n 1 ./rowmv $N
        done
done