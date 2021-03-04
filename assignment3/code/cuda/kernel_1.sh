#!/bin/bash
#SBATCH --job-name=kernel1
#SBATCH --nodes=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=30:00

module purge
module load gcc
module load cuda

nvcc -O2 nlm_kernel.cu -lm -o v1 
./v1 ../data/flower_256_noise.csv 256 1 1

