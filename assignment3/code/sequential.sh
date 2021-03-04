#!/bin/bash
#SBATCH --job-name=v_seq
#SBATCH --time=00:30:00
#SBATCH --partition=batch
#SBATCH --qos=small
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module purge
module load gcc

gcc -O2 nlm_sequencial.c -lm
 
./a.out data/flower_256_noise.csv 256 0.1 7 1.6
 
