
#include <iostream>    
#include <cmath>        
#include <algorithm>    
#include <fstream>
#include <string>
#include <iomanip>
#include <chrono>
#include <armadillo>

using namespace std;
using namespace arma;
ofstream ofile;

//Function used
double f(double);

int main()
{
    //Get desired max size of matrix A, and filename for storing execution time from user.
    
    //Gets the value n, for the max size of matrix A, from user.
    int exponent;
    cout << "Please enter the max  power 10^n, for a matrix A of size 10^n x 10^n: ";
    cin >> exponent;

    //Gets filename from user
    string filename;
    cout << "Please enter name of file you wish to write results to: ";
    cin >> filename;

    //Array for storing execution times for different sized matrices.
    double *executiontime = new double[pow(10.0, exponent)];

    //Loops the algorithm over powers of 10.
    for (int i = 1; i <= exponent; i++) {
        int n = pow(10.0, i);

        //Makes a matrix of size nxn, fills it with 0.
        mat A = zeros<mat>(n, n);
        //Updates the value of the diagonal of A to 2.
        for (uword i = 0; i < n; i++) {
            A(i, i) = 2.0;
        }
        //Updates the values above and below the diagonal of A to -1
        for (uword i = 0; i < n - 1; i++) {
            A(i, i + 1) = A(i + 1, i) = -1.0;
        }
        
        //Makes a vector b, and fills it with the values of f_i(x_i) for every i.
        vec b = vec(n);
        double h = 1.0 / n;
        for (uword i = 0; i < n; i++) {
            b(i) = (h*h * f(i * h));
        }
        
        //Fetches the time right before the algorithm starts.
        //To be used in determining the execution time.
        auto start = chrono::high_resolution_clock::now();

        //LU-decompose A, with the use of Armadillo-library. 
        mat L, U;
        lu(L, U, A);

        //Solves Ly=b, then Uv = y.
        vec y = solve(L, b);
        vec v = solve(U, y);

        //Fetches current time when algorithm is finished. 
        //Calculates the execution time, and stores it in the array Executiontime.
        auto finish = chrono::high_resolution_clock::now();
        auto timeused = chrono::duration_cast<chrono::microseconds>(finish - start).count();
        executiontime[i - 1] = timeused;
    }

    //Appends .txt to chosen filename, as to save the file as .txt-file.  
    string fileout = filename;
    fileout.append(".txt"); 
    
    //Open file with chosen name, and writes in different sized matrices that was solved used LU-decomp
    //as well as the corresponding execution time. Closes file when finished.
    ofile.open(fileout);
    ofile << "Execution times for different sized matrices [microseconds]:" << endl;
    for (int i = 0; i < exponent; i++) {
        ofile << (pow(10, i + 1)) << "x" << (pow(10, i + 1));
        ofile << setw(10) << setprecision(5) << executiontime[i] << endl;
    }
    ofile.close();
    
    return 0;
}

//Function to calculate values of f_i(x_i).
double f(double x)
{
    return 100 * exp(-10.0 * x);
}


