#!/bin/bash
#SBATCH --job-name=kernel2
#SBATCH --nodes=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=30:00

module purge
module load gcc
module load cuda


nvcc -O2 nlm_kernel_sm.cu -lm -o v2
./v2 ../data/flower_256_noise.csv 256 0.1 1.6
