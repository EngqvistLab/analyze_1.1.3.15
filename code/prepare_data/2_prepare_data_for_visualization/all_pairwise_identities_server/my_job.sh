#!/usr/bin/env bash
#SBATCH -A C3SE2018-1-24
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J all_ident
#SBATCH -t 167:00:00





# load modules
module load iccifort/2017.1.132-GCC-6.3.0-2.27 impi/2017.1.132
module load iccifort/2017.1.132-GCC-6.3.0-2.27 Python/3.6.1


# set path to find muscle executable
PATH=$PATH:~/software/muscle/

# run the job
python3 pairwise_identity_of_all.py


