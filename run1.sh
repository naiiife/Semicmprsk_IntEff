#!/bin/bash
#SBATCH --job-name=int_frailty
#SBATCH --time=12:00:00
#SBATCH --mail-type=end,fail
#SBATCH --mem=4g
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1000

R CMD BATCH --no-save --no-restore simulation1.R