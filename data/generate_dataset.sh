#!/bin/bash
#SBATCH --time 2:00:00
#SBATCH --no-requeue
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH -p debug

module load python/2.7
python generate_dataset.py
