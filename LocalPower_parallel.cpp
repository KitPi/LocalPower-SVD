#include "core/graph.h"
#include "core/utils.h"
#include <Eigen/Dense>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <atomic>
#include <random>
#include <fstream> 
#include <sstream> 
#include <vector>


uint DEFAULT_POWER_ITER =4;

uint DEFAULT_RANK =5;

using namespace std;
using namespace Eigen;

MatrixXd GS_Ortho(const MatrixXd& A) {
    //cout<< "Orthogonalising ..." << endl;
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
    //cout<<"done"<<endl;
    return Q;
}


// Function to read a CSV file into an Eigen matrix 
MatrixXd readCSV(const string& file) { 
    ifstream infile(file); 
    if (!infile.is_open()) { 
        throw runtime_error("Could not open file"); 
    } 

    vector<vector<double>> data;
    string line; 
    // Skip the first line (headers) 
    getline(infile, line);
    //cout << line<<endl;
    while (getline(infile, line)) { 
        istringstream iss(line); 
        string value;
        vector<double> row; 
        
        // Read each value, separated by commas (,)
        while (getline(iss, value, ',')) { 
            row.push_back(stod(value)); 
        } 
        data.push_back(row); 
    } 
    infile.close(); 
    // Determine the dimensions of the matrix 
    
    int rows = data.size(); 
    int cols = data.empty() ? 0 : data[0].size(); 
    //cout << rows <<", "<< cols << endl;
    // Creeate an Eigen matrix and populate it with the values 
    MatrixXd matrix(rows, cols); 
    for (int i = 0; i < rows; ++i) { 
        for (int j = 0; j < cols; ++j) { 
            matrix(i, j) = data[i][j]; 
        }
    } 
    return matrix;
}

// Barriers 
CustomBarrier* barrier1;
CustomBarrier* barrier2;

//Atomics


void LocalPower(MatrixXd* A, MatrixXd* Z_1, atomic<double>* Y, double* time, uint p, uint T, uint k, uint nThread, uint start, uint end) {
    // Get A_i and A_i transposed
    MatrixXd A_i =A->block(start, 0, (end - start), A->cols());
    MatrixXd A_iT = A_i.transpose();
    // Get Y_i 
    //cout<< "making Y_i" << endl;
    MatrixXd Y_i = MatrixXd::Zero(A->cols(), k);
    for (int i = 0; i< A->cols(); i++){
        for (int j = 0; j< k; j++){
            Y_i(i, j) = Y[i*k + j].load();
        }
    }

    // Get D_i
    MatrixXd D_i = MatrixXd::Identity(k,k);

    timer t1;
    t1.start();

    for (int t=0; t<T; t++){

        // Perform GS orthogonalisation 
        MatrixXd Z_i = GS_Ortho(Y_i);

        if(((t % p) == p-1) || (t==T-1)){
            // Timer increment
            *time += t1.stop();

            // Final iteration condition: Calculate D_i and Y_i, then aggregate add Y
            if(nThread == 0){
                // Communication
                *Z_1=Z_i;
            }

            //Barrier: Wait for thread 0 to send Z_1
            barrier1->wait();

            if(nThread>0){
                // Find D_i for threads > 0
                // (Sign fixing): D_i[j,j] = sign(<Z_i[:,j],Z_1[:,j]>)
                for (int j=0; j < k; j++){
                    D_i(j,j) = (Z_i.col(j).dot(Z_1->col(j)) >= 0) ? 1.0: -1.0;
                }
            }
            // (Final iteration of Y_i): Y_i = (1/s) *A_i^T *A_i *Z_i *D_i
            Y_i = 1.0/(double)(end - start) * (A_iT * (A_i)) * (Z_i) * (D_i);
            
            // -- Linearization Point
            // atomic add : Y += p * Y_i 
            double expected;
            double new_val;
            bool swap;
            for (int i =0; i< A->cols(); ++i){
                for(int j=0; j< k; ++j){
                    swap = false;
                    expected = Y[i*k +j].load();
                    double frac = (double)(end-start)/(double)A->rows();
                    new_val = expected+ frac *Y_i(i, j);
                    // Atomic compare and swap weak: ordering doesn't matter because Y is sum-accumulating.
                    do{
                        swap = Y[i*k +j].compare_exchange_weak(expected, new_val);
                        if(!swap) {
                            new_val = expected+ frac *Y_i(i, j);
                        }
                    } while(!swap);
                }
            }

            // -- Linearization Point
            //Barrier: Homogeneity of Y
            barrier2->wait();
            // Get new Y_i 
            for (int i = 0; i< A->cols(); i++){
                for (int j = 0; j< k; j++){
                    Y_i(i, j) = Y[i*k + j].load();
                }
            }
            // Timer start again
            t1.start();
        }else{
            //Non Aggregative iteration O(kds)
            Y_i = 1.0/(double)(end - start) * (A_iT * A_i) * (Z_i);
        }
    }
}

int main (int argc, char *argv[]) {
    cxxopts::Options options(
        "Local Power",
        "Calculate the k-th Truncated Singular Vectors of the input matrix, in parallel");
    options.add_options(
        "",
        {
            {"nThreads", "Number of Threads",
             cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
            {"nIterations", "Maximum number of iterations",
             cxxopts::value<uint>()->default_value(to_string(32))},
            {"powerIter", "Number of power iterations",
             cxxopts::value<uint>()->default_value(to_string(DEFAULT_POWER_ITER))},
            {"rank", "Rank of the truncated approximate",
             cxxopts::value<uint>()->default_value(to_string(DEFAULT_RANK))},
            {"inputFile", "Input csv file path",
             cxxopts::value<string>()->default_value(
               "/home/kpeer/CMPT_431/LocalPower-SVD/data/covtype.csv")},
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
    // Input csv file
    string input_file_path = cl_options["inputFile"].as<string>();

    cout << fixed;
    cout << "Number of Threads : " << n_threads << endl;
    cout << "Total iterations : " << T << endl;
    cout << "Power iterations: " << p << endl;
    cout << "Rank: " << k << endl;
    cout << "CSV file: " << input_file_path << endl;


    // Create barriers
    barrier1= new CustomBarrier(n_threads);
    barrier2= new CustomBarrier(n_threads);

    // start global timer
    double time_taken = 0.0;
    timer t1;
    t1.start();

    //Import Matrix A from .csvf
    MatrixXd A= readCSV(input_file_path);

    // Start a random distribution for inititalising Y
    // Create a random number generator 
    random_device rd; 
    mt19937 gen(rd()); 
    normal_distribution<> d(0.0, 1.0); // mean 0, standard deviation 1


    // Instantiate Y as a gaussian matrix (cols,k)
    double value;
    int size = A.cols()*k;
    std::vector<std::atomic<double>> Y(size);// = {0.0};
    for (int i =0; i< A.cols(); i++){
        for (int j =0; j<k; j++){
            value = d(gen);
            Y[i*k + j].store((double)value);
        }
    }
    
    // Times
    double times[n_threads] ={0.0};

    //cout<< "Generating Start points"<<endl;
    // Split the Matrix 
    int start_points[n_threads+1];  
    int start=0;
    for (int i=0; i< n_threads+1;i++){
        start_points[i]=start;
        start+=(i < A.rows() % n_threads) ? int(A.rows() / n_threads) + 1 : int(A.rows()/n_threads);
    } 

    // Z_1
    MatrixXd Z_1(A.cols(), k);

    // Start threads
    //cout<< "Starting Threads: "<<endl;
    thread threads[n_threads];
    for (int i=0; i< n_threads; ++i){
        threads[i]=thread(LocalPower, (MatrixXd*) &A, (MatrixXd*) &Z_1, (atomic<double>*) Y.data(), (double*) &times[i], (uint) p, (uint) T, (uint) k,(uint) i, (uint) start_points[i], (uint) start_points[i+1]);
        //MatrixXd* A, MatrixXd* Z_1, atomic<double>* Y, double* time, uint p, uint T, uint k, uint nThread, uint start, uint end)
    }

    // Join threads
    for (int i=0; i<n_threads;++i){
        threads[i].join();
    }
    time_taken = t1.stop();

    cout << "Cumulative computation time between communication intervals:" << endl;
    for (int i=0;i<n_threads;i++){
        cout << i << ", "<< times[i] << endl;
    }

    cout << "Total time including communication intervals:" << endl;
    cout << time_taken <<endl;

    //cout<<A.block(0,0,k,k)<<endl;
    MatrixXd Yf = MatrixXd::Zero(A.cols(), k);

    for(int i =0; i<A.cols(); i++){
        for(int j=0; j<k; j++){
            Yf(i, j)=Y[i*k+j].load();
        }
    }

    
    MatrixXd Z = GS_Ortho(Yf);

    MatrixXd err = readCSV("data/matrix.csv");
    MatrixXd diff = Z - err;
    cout << "Error: " << diff.norm() <<endl<<endl;

    //cout<< Z.rows() << ", "<< Z.cols() <<endl;
    //cout<< Z<<endl;
    return 0;

}