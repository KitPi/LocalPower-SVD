#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --time=5:00
#SBATCH --partition=fast


#SBATCH --ntasks=4
srun ./LocalPower_distributed --nIterations 16 --powerIter 8 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 8 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 8 --rank 8
#
#srun ./LocalPower_distributed --nIterations 16 --powerIter 4 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 4 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 4 --rank 8
#
#srun ./LocalPower_distributed --nIterations 16 --powerIter 2 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 2 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 2 --rank 8


##SBATCH --ntasks=4
#srun ./LocalPower_distributed --nIterations 16 --powerIter 8 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 8 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 8 --rank 8
#
#srun ./LocalPower_distributed --nIterations 16 --powerIter 4 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 4 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 4 --rank 8
#
#srun ./LocalPower_distributed --nIterations 16 --powerIter 2 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 2 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 2 --rank 8
#
#
##SBATCH --ntasks=2
#srun ./LocalPower_distributed --nIterations 16 --powerIter 8 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 8 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 8 --rank 8
#
#srun ./LocalPower_distributed --nIterations 16 --powerIter 4 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 4 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 4 --rank 8
#
#srun ./LocalPower_distributed --nIterations 16 --powerIter 2 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 2 --rank 8
#srun ./LocalPower_distributed --nIterations 16 --powerIter 2 --rank 8
