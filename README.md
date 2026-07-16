To run the serial SVD, please use
```Bash
./LocalPower_serial 8
```
To run the thread parallel version of LocalPower, please use
```Bash
./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 8 --rank 8
```
To run the distributed version of LocalPower, please use
```Bash
sbatch slurm/LP_dist.sh
```
