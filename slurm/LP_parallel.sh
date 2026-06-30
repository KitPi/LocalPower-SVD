#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --time=10:00
#SBATCH --partition=fast

#SBATCH --ntasks=1


srun ./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 8 --rank 8 
srun ./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 8 --rank 8 
srun ./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 8 --rank 8 

srun ./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 4 --rank 8 
srun ./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 4 --rank 8 
srun ./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 4 --rank 8 

srun ./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 2 --rank 8 
srun ./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 2 --rank 8 
srun ./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 2 --rank 8 


srun ./LocalPower_parallel --nThreads 4 --nIterations 16 --powerIter 8 --rank 8 
srun ./LocalPower_parallel --nThreads 4 --nIterations 16 --powerIter 8 --rank 8 
srun ./LocalPower_parallel --nThreads 4 --nIterations 16 --powerIter 8 --rank 8 

srun ./LocalPower_parallel --nThreads 4 --nIterations 16 --powerIter 4 --rank 8 
srun ./LocalPower_parallel --nThreads 4 --nIterations 16 --powerIter 4 --rank 8 
srun ./LocalPower_parallel --nThreads 4 --nIterations 16 --powerIter 4 --rank 8 

srun ./LocalPower_parallel --nThreads 4 --nIterations 16 --powerIter 2 --rank 8 
srun ./LocalPower_parallel --nThreads 4 --nIterations 16 --powerIter 2 --rank 8 
srun ./LocalPower_parallel --nThreads 4 --nIterations 16 --powerIter 2 --rank 8 


srun ./LocalPower_parallel --nThreads 2 --nIterations 16 --powerIter 8 --rank 8 
srun ./LocalPower_parallel --nThreads 2 --nIterations 16 --powerIter 8 --rank 8 
srun ./LocalPower_parallel --nThreads 2 --nIterations 16 --powerIter 8 --rank 8 

srun ./LocalPower_parallel --nThreads 2 --nIterations 16 --powerIter 4 --rank 8 
srun ./LocalPower_parallel --nThreads 2 --nIterations 16 --powerIter 4 --rank 8 
srun ./LocalPower_parallel --nThreads 2 --nIterations 16 --powerIter 4 --rank 8 

srun ./LocalPower_parallel --nThreads 2 --nIterations 16 --powerIter 2 --rank 8 
srun ./LocalPower_parallel --nThreads 2 --nIterations 16 --powerIter 2 --rank 8 
srun ./LocalPower_parallel --nThreads 2 --nIterations 16 --powerIter 2 --rank 8 
