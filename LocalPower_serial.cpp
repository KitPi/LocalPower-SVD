#include <Eigen/Dense>
#include <iostream>
#include "core/utils.h"

MatrixXd GS_Ortho(const MatrixXd& A) {
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

int main () {
    // Initialise matrix A
    MatrixXd A;

    // Number of singular values to approximate
    int k = 5;

    // start global timer
    timer t1;
    t1.start();

    // Compute the SVD decomposition
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);

    // Compute the truncated SVD matrix: A_k: O(n^2 k)
    MatrixXd U_k = svd.matrixU().leftCols(k);
    MatrixXd S_k = svd.singularValues().head(k).asDiagonal();
    MatrixXd V_kT = svd.MatrixV().leftCols(k).transpose();
    MatrixXd A_k = U_k * S_k * V_kT;
}
