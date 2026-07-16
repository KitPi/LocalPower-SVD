#include "core/graph.h"
#include "core/utils.h"
#include <Eigen/Dense>
#include <iostream>
#include <stdlib.h>
#include <random>
#include <fstream> 
#include <sstream> 
#include <vector>
#include <mpi.h>

uint DEFAULT_POWER_ITER =4;
uint DEFAULT_RANK =5;

using namespace Eigen;
using namespace std;

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

MatrixXd LocalPower_distributed(MatrixXd A, MatrixXd Y, double frac, uint p, uint T, uint k, MPI_Comm comm){
    //A_i, Y, frac, p, T, k, comm
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    //Matrix M
    //cout << "Matrix m: " << world_rank << endl;
    MatrixXd M = (1.0/(double)A.rows()) * A.transpose() * A;
    //if(world_rank == 0) cout<< M <<endl;

    MatrixXd At = A.transpose();

    //cout << "Matrix Z: " << world_rank << endl;
    MatrixXd Z_1(A.cols(), k);
    MatrixXd Z(A.cols(), k);

    // Get D_i
    //cout << "Matrix D: " << world_rank << endl;
    MatrixXd D_i = MatrixXd::Identity(k,k);

    //cout<< frac << ", " << world_rank<<endl;

    //Timer
    double time = 0.0;
    timer t1;
    t1.start();

    for (int t=0; t<T; t++) {
        // Perform GS orthogonalisation 

        //cout<< "Iteration: " << t << ", " << world_rank <<endl;

        //cout << "Orthogonalise: " << world_rank << endl;
        Z = GS_Ortho(Y);

        if(((t % p) == p-1) || (t==T-1)){
            //if (world_rank == 0) cout<< Y << endl;
            //Iterate time 
            time += t1.stop();

            // Communication tsarts here
            if(world_rank == 0){
                // Send Z_1 to all other processes
                Z_1 = Z;
                // Braodcast Z_1 from process 0 to all other processes
                //cout<< "Broadcast from process 0"<<endl;
                //MPI_Bcast(Z_1.data(), A.cols() * k, MPI_DOUBLE, 0, comm);
            } 

            //Barrier
            //cout<< "Barrier 1: " << world_rank << endl;
            MPI_Barrier(comm);

            MPI_Bcast(Z_1.data(), A.cols() * k, MPI_DOUBLE, 0, comm);
            
            for (int j=0; j < k; j++){
                D_i(j,j) = (Z.col(j).dot(Z_1.col(j)) >= 0) ? 1.0: -1.0;
            }
            
            
            Y = frac * M * (Z) * (D_i);
            //if(world_rank == 0) cout<< Y <<endl;

            //vector<double> globalY;
            //globalY.resize(A.cols() * k);
            //MPI_Allreduce(localY.data(), globalY.data(), A.cols()*k, MPI_DOUBLE, MPI_SUM, comm)
            MatrixXd buffer(A.cols(), k);
            MPI_Allreduce(Y.data(), buffer.data(), A.cols()*k, MPI_DOUBLE, MPI_SUM, comm);
            Y= buffer;
            

            // Barrier
            //cout<< "Barrier 2: " << world_rank << endl;
            MPI_Barrier(comm);

            //timer
            t1.start();
            //cout<< Y.block(0,0,k,k)<< " "<<world_rank << endl;

            

        } else {
            Y = 1.0/(double)(A.rows()) * (At * A) * Z;
        }
    }

    vector<double> time_array;
    time_array.resize(world_size);
    MPI_Gather(&time, 1, MPI_DOUBLE, time_array.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (world_rank == 0) {
        cout<< "Cumulative computation time between communication intervals:"<<endl;
        for (int i=0; i<time_array.size(); i++) {
            cout<< i << ", " << time_array[i] << endl;
        }
    }
    //MPI_Allreduce(localY, globalY, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD)
    //MatrixXd* Z_1,
    return GS_Ortho(Y);
}

int main (int argc, char *argv[]) {
    //Initiatlise processes
    MPI_Init(NULL, NULL);
    
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    //Command line args
    cxxopts::Options options(
        "Local Power",
        "Calculate the k-th Truncated Singular Vectors of the input matrix, distributed");
    options.add_options(
        "",
        {
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
    // Number of iterations
    uint T = cl_options["nIterations"].as<uint>();
    // Number of power iterations
    uint p = cl_options["powerIter"].as<uint>();
    // Number of singular values to approximate
    uint k = cl_options["rank"].as<uint>();
    // Input csv file
    string input_file_path = cl_options["inputFile"].as<string>();

    if(world_rank == 0){
        cout << fixed;
        cout << "Total iterations : " << T << endl;
        cout << "Power iterations: " << p << endl;
        cout << "Rank: " << k << endl;
        cout << "CSV file: " << input_file_path << endl;
        cout << "World size: " << world_size << endl;
    }
    

    // Initialise matrix A
    MatrixXd A = readCSV(input_file_path);

    // Start a random distribution for inititalising Y
    // Create a random number generator 
    random_device rd; 
    mt19937 gen(rd()); 
    normal_distribution<> d(0.0, 1.0); // mean 0, standard deviation 1

    // Instantiate Y as a gaussian matrix (cols,k)
    //double value;
    //int size = A.cols()*k;
    //double Y[size];// = {0.0};
    double value;
    MatrixXd Y(A.cols(), k);
    for (int i =0; i< A.cols(); i++){
        for (int j =0; j<k; j++){
            value = d(gen);
            Y(i, j) = (double)value;
        }
    }

    // Times
    double times[world_size] ={0.0};

    // start global timer
    double time_taken = 0.0;
    timer t1;
    t1.start();

    // Send A submatrices to other processes:
    MatrixXd A_i;
    double frac;
    if (world_rank == 0){
        //cout<< "Loading Matrix from CSV: (process 0)"<<endl;
        MatrixXd A = readCSV(input_file_path);
        int numRows = A.rows() / world_size;

        // Generating starting points
        //cout<< "Generating starting points" <<endl;
        int start_points[world_size+1];  
        int start=0;
        for (int i=0; i< world_size+1;i++){
            start_points[i]=start;
            start+=(i < A.rows() % world_size) ? int(A.rows() / world_size) + 1 : int(A.rows()/world_size);
        }

        //Send A matrices
        //cout<< "Send A matrices" << endl;
        for (int i = 1; i < world_size; i++) {
            MatrixXd m = A.block(start_points[i], 0, start_points[i+1] - start_points[i], A.cols());
            int rows = m.rows();
            int cols = m.cols();
            frac = (double)(m.rows()) / (double)A.rows(); 
            MPI_Send(&rows, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&cols, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(m.data(), rows * cols, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
            MPI_Send(&frac, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
        }
        A_i = A.block(start_points[0], 0, start_points[1], A.cols());
        frac = (double)(start_points[1])/(double)A.rows();
        //cout << "Matrices sent" << endl;
    } else {
        // Receive Matrix A_i from processs 0;
        int rows, cols;
        MPI_Recv(&rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&cols, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        A_i.resize(rows, cols);
        MPI_Recv(A_i.data(), rows * cols, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&frac, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //cout << "Receive A Matrix: " << world_rank << endl;
    }

    // LocalPower
    //cout << "starting: " << world_rank << endl;
    MatrixXd r = LocalPower_distributed(A_i, Y, frac, p, T, k, MPI_COMM_WORLD);
    //cout << "ending: " << world_rank << endl;
    // MatrixXd A, uint p, uint T, uint k, MPI_Comm comm

    if (world_rank == 0) {
        cout<< "Total time including communication intervals:" <<endl;
        cout<< t1.stop() <<endl;

        MatrixXd err = readCSV("data/matrix.csv");
        MatrixXd diff = r - err;
        cout << "Error: " << diff.norm() <<endl<<endl;
        //cout<< r << endl;
    }
    // Finaslise the distributed processes
    MPI_Finalize();

    return 0;

}