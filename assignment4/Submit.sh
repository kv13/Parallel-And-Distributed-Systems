#!/bin/bash
#SBATCH --job-name=d_s_2_2
#SBATCH --time=00:10:00
#SBATCH --partition=batch
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1

echo shared memory and distributed memory optimization

module purge
module load gcc
module load openmpi

mpicc -c src/mmio.c         -O3 -o src/mmio.o
mpicc -c src/utils.c        -O3 -o src/utils.o 
mpicc -c src/blocking.c     -O3 -o src/blocking.o
mpicc -c src/bmm.c          -O3 -o src/bmm.o
mpicc -c src/parallel_bmm.c -O3 -o src/parallel_bmm.o -fopenmp

mpicc src/mmio.o src/utils.o src/blocking.o src/bmm.o src/parallel_bmm.o blocked_bmm_distmem.c -O3 -fopenmp -o script_4

export OMP_NUM_THREADS=$SLURM_NTASKS

echo bmm
srun -n 2 ./script_4 data/matrix_A_100k.mtx data/matrix_B_100k.mtx 10000
srun -n 2 ./script_4 data/matrix_A_300k.mtx data/matrix_B_300k.mtx 75000

echo bmm filtered
srun -n 2 ./script_4 data/matrix_F_1M.mtx data/matrix_A_1M.mtx data/matrix_B_1M.mtx 250000
srun -n 2 ./script_4 data/matrix_F_3M.mtx data/matrix_A_3M.mtx data/matrix_B_3M.mtx 375000
srun -n 2 ./script_4 data/matrix_F_5M.mtx data/matrix_A_5M.mtx data/matrix_B_5M.mtx 250000
