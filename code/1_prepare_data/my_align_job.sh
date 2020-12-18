#!/usr/bin/env bash
#SBATCH -A C3SE2020-1-14
#SBATCH -n 64
#SBATCH -J align
#SBATCH -t 0-12:00 # time (D-HH:MM)

module load GCC/8.3.0  OpenMPI/3.1.4
module load Biopython/1.75-Python-3.7.4
export PATH=$PATH:/cephyr/users/marengq/Vera/

python align_for_server.py $SLURM_ARRAY_TASK_ID
echo 'script finished'
