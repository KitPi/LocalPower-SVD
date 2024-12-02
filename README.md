# LocalPower-SVD
Implementation and analysis of of methods for distributed/parallel Singular Value Decomposition using local power iteration. 
distributed/parallel SVD [[2002.08014v4.pdf]]

# Todo
- [ ] [[CSIL labs]]
- [ ] https://coursys.sfu.ca/2021fa-cmpt-470-d1/pages/Docker
- [ ] https://ieeexplore.ieee.org/document/7868429 https://github.com/NLKNguyen/alpine-mpich/wiki/Multi-Host-Orchestration

## Computation Time
 **(M workers vs. P iterations) vs. Computation Time**
 Reduced communication efficiency with smaller $P$ iterations
 $O(dkm)$ communication complexity : for each communication, after P iterations
## Accuracy 
 **(M Workers vs. P iterations) vs. Error**
 Reduced accuracy with larger values of $P$ and for smaller granularity of sub-matrices (More Workers $M$)

## Sign fixing on Local Step
 Reduces the complexity of local workers, 
 $D^{(i)}_t [j, j] = sgn ( 〈 Z^{(i)} _t [:, j] , Z^{(1)} _t [:, j] 〉 ) , ∀ j ∈ [k]$
