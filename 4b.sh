#!/bin/bash
# Job name:
#SBATCH --job-name=aspen_ipm_4b
# Account:
#SBATCH --account=fc_mel
#
# Partition:
#SBATCH --partition=savio3
#
# Request one node:
#SBATCH --nodes=1
#
# Specify one task:
#SBATCH --ntasks-per-node=1
#
# Number of processors for threading:
#SBATCH --cpus-per-task=40
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
## Command(s) to run (example):
module load r
module load r-spatial
R CMD BATCH --no-save 4b.R aspen_ipm_4b.Rout