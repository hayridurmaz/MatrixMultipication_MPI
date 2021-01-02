#!/bin/bash
###################################################################################################
#SBATCH --job-name=1_hayridurmaz-rowmv
#SBATCH --account=users
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=32000
#SBATCH --partition=cpu
#SBATCH --output=1_hayridurmaz-rowmv.txt
###################################################################################################

module load gcc/10.2.0
#echo $HOSTNAME
#echo $SLURM_CPUS_PER_TASK
make
#cat /proc/cpuinfo
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export NB=256
for N in 2048 4096 8192 16384; do
        for nthreads in 1 2 4 8 16 32; do
                export OMP_NUM_THREADS=$nthreads
                srun  -n 1 ./rowmv -n $N -nb $NB
        done
done