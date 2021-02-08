#!/usr/bin/env bash
# invoke using sbatch my_filter_job.sh
#SBATCH -A C3SE2020-1-14
#SBATCH -n 64
#SBATCH -J filter
#SBATCH -t 6-12:00 # time (D-HH:MM)

module load GCC/8.3.0  OpenMPI/3.1.4
#module load Biopython/1.75-Python-3.7.4
module load Anaconda3/5.3.0

# activate the pre-installed environment
source activate hmmer

export PATH=$PATH:/cephyr/users/marengq/Vera/

python3 filter_for_server.py
echo 'script finished'
