#!/usr/bin/env bash
# invoke using sbatch --array=1-7 my_hmm_job.sh
#SBATCH -A C3SE2020-1-14
#SBATCH -n 64
#SBATCH -N 1
#SBATCH -J hmmer
#SBATCH -t 2-12:00 # time (D-HH:MM)

module load GCC/8.3.0  OpenMPI/3.1.4
#module load Biopython/1.75-Python-3.7.4
module load Anaconda3/5.3.0

# activate the pre-installed environment
source activate hmmer

export PATH=$PATH:/cephyr/users/marengq/Vera/

echo $SLURM_ARRAY_TASK_ID
python3 run_hmmsearch.py $SLURM_ARRAY_TASK_ID
echo 'script finished'
