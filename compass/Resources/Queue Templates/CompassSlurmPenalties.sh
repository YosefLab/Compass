#!/bin/sh
#SBATCH --job-name=compass_penalties # Job name
#SBATCH --ntasks={MAX_PROCESSES}    # Run a single task	
#SBATCH --cpus-per-task=1           # Number of CPU cores per task
#SBATCH --output=array_%A-%a.log    # Standard output and error log

CONFIG=$1

# --data and other options is taken from the config file
compass --only-penalties --config-file $CONFIG
