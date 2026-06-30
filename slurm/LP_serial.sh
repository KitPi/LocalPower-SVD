#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --time=2:00
#SBATCH --partition=fast



#SBATCH --ntasks=1
srun ./LocalPower_serial 8
