#!/usr/bin/env bash
#SBATCH -A C3SE2020-1-14
#SBATCH -n 64
#SBATCH -J ktup
#SBATCH -t 2-12:00 # time (D-HH:MM)
module load GCC/8.3.0  OpenMPI/3.1.4
module load Biopython/1.75-Python-3.7.4
#pip install --user alfpy biopython
echo 'script starting'
python ktup_for_server.py $SLURM_ARRAY_TASK_ID
echo 'script finished'
