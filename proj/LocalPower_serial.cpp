#include "core/utils.h"
#include <Eigen/Dense>
#include <iostream>
#include <stdlib.h>
//#include <cstdlib>
#include <fstream> 
#include <sstream> 
#include <vector>
using namespace std;
using namespace Eigen;

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

int main (int argc, char* argv[]) {
    // Initialise matrix A
    MatrixXd A = readCSV("data/covtype.csv");
    MatrixXd err = readCSV("data/matrix.csv");

    int k = stoi(argv[1]);
    cout << "Rank: " << k << endl;

    // Number of singular values to approximate
    //int k = (8);

    // start global timer
    timer t1;
    t1.start();

    // Compute the SVD decomposition
    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);

    // Compute the truncated SVD matrix: A_k: O(n^2 k)
    //MatrixXd U_k = svd.matrixU().leftCols(k);
    //MatrixXd S_k = svd.singularValues().head(k).asDiagonal();
    MatrixXd V_k = svd.matrixV().leftCols(k);
    //MatrixXd A_k = U_k * S_k * V_k.transpose();

    double time = t1.stop();

    MatrixXd diff = V_k - err;
    cout << "Total time:"<< time <<endl;
    cout << "Error: " << diff.norm() <<endl;

    //cout<< V_k.rows() << ", "<< V_k.cols() <<endl;
    //cout<< V_k <<endl;
}