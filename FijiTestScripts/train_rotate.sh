#!/bin/bash
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --partition=nvidia-a100
#SBATCH --time=168:00:00
#SBATCH --ntasks=16
#SBATCH --output=./train_r_19-%j-%N.out
#SBATCH --error=./train_r_19-%j-%N.err
#SBATCH --job-name=train_r_19
#SBATCH --gres=gpu:1

python Scripts/create_optimized_rotate.py
