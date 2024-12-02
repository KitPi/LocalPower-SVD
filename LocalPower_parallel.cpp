#include "core/graph.h"
#include "core/utils.h"
#include <Eigen/Dense>
#include <iostream>
#include <thread>

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

// Barrier 
CustomBarrier barrier1;
CustomBarrier barrier2;

void LocalPower(MatrixXd* A, MatriXd* Y, uint p, uint T, uint start, uint end) {
    MatrixXd A_i =A->cols(start, end);
    MatrixXd Y_i =Y->cols(start, end);

    for (int i=0; i<T, i++){
        MatrixXd Z_i = 
    }
}

int main (int argc, char *argv[]) {
    cxxopts::Options options(
        "Local Power",
        "Calculate the Singular Value Decomposition of the input matrix, in parallel");
    options.add_options(
        "",
        {
            {"nThreads", "Number of Threads",
             cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
            {"nIterations", "Maximum number of iterations",
             cxxopts::value<uint>()->default_value(32)},
             {"powerIter", "Number of power iterations",
             cxxopts::value<uint>()->default_value(8)},
            {"rank", "Rank of the truncated approximate",
             cxxopts::value<uint>()->default_value(5)},
        });

    auto cl_options = options.parse(argc, argv);
    // Number of threads
    uint n_threads = cl_options["nThreads"].as<uint>();
    // Number of iterations
    uint T = cl_options["nIterations"].as<uint>();
    // Number of power iterations
    uint p = cl_options["powerIter"].as<uint>();
    // Number of singular values to approximate
    uint k = cl_options["rank"].as<uint>();

    // Create barriers
    barrier1= new CustomBarrier(n_threads);
    barrier2= new CustomBarrier(n_threads);

    // start global timer
    timer t1;
    t1.start();

    // Start threads
    std::thread threads[n_threads];
    for (int i=0; i< n_threads; ++i){
      threads[i]=std::thread();
    }

    // Join threads
    for (int i=0; i<n_threads;++i){
      threads[i].join();
    }
    time_taken = t1.stop();

}