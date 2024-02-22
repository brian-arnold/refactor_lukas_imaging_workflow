#!/bin/bash
#SBATCH -J template
#SBATCH -o out_part1
#SBATCH -e err_part1
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=22        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=6G         # memory per cpu-core (4G is default)
#SBATCH --time 1-00:00:00        # DAYS-HOURS:MINUTES:SECONDS

source ~/miniforge3/etc/profile.d/conda.sh
conda activate ZZ_calcium


python analysis_pipeline_orco_part1.py
