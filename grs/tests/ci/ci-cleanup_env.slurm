#!/bin/bash
#
#SBATCH --job-name=grs-cleanup_env
#SBATCH --time=00:30:00
#SBATCH -N 1
#SBATCH --account=campus

#cnes_level2

ml conda
conda env remove -p /work/scratch/data/${USER}/.conda/envs/grs2_tests
