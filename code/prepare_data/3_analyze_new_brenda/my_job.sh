#!/usr/bin/env bash
#SBATCH -A C3SE2020-1-14
#SBATCH -n 64
#SBATCH -J ktup
#SBATCH -t 3-00:00 # time (D-HH:MM)

module load GCC/8.2.0-2.31.1  CUDA/10.1.105
module load OpenMPI/3.1.3

python3 ktup_for_server.py
echo 'script finished'
