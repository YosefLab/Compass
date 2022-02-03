#!/bin/sh
#SBATCH --job-name=compass_array # Job name
#SBATCH --ntasks=1                   # Run a single task	
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --output=array_%A-%a.log    # Standard output and error log
#SBATCH --array=1-{NUM_SAMPLES}%{MAX_PROCESSES}          # Array range

CONFIG=$1

# --data and other options is taken from the config file
compass --single-sample $SLURM_ARRAY_TASK_ID --config-file $CONFIG
