#!/bin/bash
#SBATCH --job-name=spatial_ccc_experiments
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
#SBATCH --time=10:00:00
#SBATCH --output=logs/main_%j.out
#SBATCH --error=logs/main_%j.out
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_normal
#SBATCH --nice=1000

python /home/icb/francesca.drummer/1-Projects/spatial_ccc_experiments/src/methods/commot/commot_script.py

