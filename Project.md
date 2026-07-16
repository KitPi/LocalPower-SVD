[[CMPT 431]]

# LocalPower-SVD
Implementation and analysis of of methods for distributed/parallel Singular Value Decomposition using local power iteration. 
distributed/parallel SVD [[2002.08014v4.pdf]]

For $t=[0,T]:$ 
	**Broadcast:** If $t\in\mathbf I_T$, the server sends $Y_t$ to workers; let $Y_t^{(i)}\leftarrow Y_{(t)};$ 
	**Local Computation:** For all $i\in[m]$, the $i^{th}$ worker locally computes $Z_t^{(i)}=orth(Y_t^{(i)})$ and $Y_{t+1}^{(i)}=\frac{1}{s_i}A_i^TA_iZ_t^{(i)};$ 
	**Aggregation:** If $(t+1)\in\mathbf I_T,$ the server computes $D^{(i)}_t [j, j] = sgn ( 〈 Z^{(i)} _t [:, j] , Z^{(1)} _t [:, j] 〉 ) , ∀ j ∈ [k];$ $Y^{(i)}_{(t+1)}= \frac{1}{s_i}A_i^TA_iZ_t^{(i)}D_t^{(i)};$ $Y_{(t+1)}=\sum\limits_{i=1}^{m}p_iY^{(i)}_t;$ 

## Computation Time and Communication Time
 **(M workers vs. P iterations) vs. Computation Time**
 Reduced communication efficiency with smaller $P$ iterations
 $O(dkm)$ communication complexity : for each communication, after P iterations
## Accuracy 
 **(M Workers vs. P iterations) vs. Error**
 Reduced accuracy with larger values of $P$ and for smaller granularity of sub-matrices (More Workers $M$)

## Sign fixing on Local Step
 Reduces the complexity of local workers, as a substitute for OPT SVD
 At iteration $t=P-1$ 
 $D^{(i)}_t [j, j] = sgn ( 〈 Z^{(i)} _t [:, j] , Z^{(1)} _t [:, j] 〉 ) , ∀ j ∈ [k]$
 $Y^{(i)}_{(t+1)}= \frac{1}{s_i}A_i^TA_iZ_t^{(i)}D_t^{(i)}$
 


## Eigen
 git clone https://gitlab.com/libeigen/eigen.git


### Orthogonal
```C++
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

MatrixXd gramSchmidt(const MatrixXd& A) {
    int m = A.rows();
    int n = A.cols();
    MatrixXd Q = MatrixXd::Zero(m, n);
    VectorXd u, q;
    for (int j = 0; j < n; ++j) {
        u = A.col(j);
        for (int i = 0; i < j; ++i) {
            q = Q.col(i);
            u -= (q.dot(A.col(j)) / q.dot(q)) * q;
        }
        Q.col(j) = u / u.norm();
    }
    return Q;
}

int main() {
    // Example matrix A
    MatrixXd A(3, 2);
    A << 1, 1,
         1, 0,
         0, 1;
    MatrixXd Q = gramSchmidt(A);
    cout << "Orthogonalized matrix Q:\n" << Q << endl;
    return 0;
}
```

### SVD
```C++
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main() {
    // Example matrix A
    MatrixXd A(3, 3);
    A << 1, 2, 3,
         4, 5, 6,
         7, 8, 9;
    
    // Compute the SVD
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    
    // Number of singular values to keep
    int k = 2;
    
    // Truncate the SVD
    MatrixXd U_k = svd.matrixU().leftCols(k);
    MatrixXd S_k = svd.singularValues().head(k).asDiagonal();
    MatrixXd V_kT = svd.matrixV().leftCols(k).transpose();
    
    // Form the approximate matrix A_k
    MatrixXd A_k = U_k * S_k * V_kT;
    
    // Output the results
    cout << "Original matrix A:\n" << A << endl;
    cout << "\nApproximated matrix A_k:\n" << A_k << endl;
    
    return 0;
}

```

```Bash
./LocalPower_serial 8
```
```Bash
./LocalPower_parallel --nThreads 8 --nIterations 16 --powerIter 8 --rank 8
```
```Bash
./LocalPower_distributed --nIterations 16 --powerIter 8 --rank 8
```


https://www.kaggle.com/datasets/uciml/forest-cover-type-dataset/data