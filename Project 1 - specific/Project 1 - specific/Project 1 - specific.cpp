
#include <iostream>    
#include <cmath>        
#include <algorithm>    
#include <fstream>
#include <string>
#include <iomanip>
#include <chrono>

using namespace std;
ofstream ofile;

//Functions used
double f(double);
double f_exact(double);

int main()
{
    
    //Fetch max power of 10 from user.
    int exponent;
    cout << "Please enter the max  power 10^n: ";
    cin >> exponent;

    //Fetch filename from user.
    string filename;
    cout << "Please enter name of file you wish to write results to.\n";
    cout << "File will be stored as a .txt-file, and the name will have the _10^n added to it, for each n: ";
    cin >> filename;

    //Arrays to store execution time and max relative error for different 10^n.
    //Used to write that to a seperate file.
    double* ExecutionTimes = new double[exponent];
    double* MaxRelativeError = new double[exponent];

    //Loops the algorithm over powers of 10.
    for (int i = 1; i <= exponent; i++) {
        int n = pow(10.0, i);

        //Attachs _10^n.txt for as to make a file for each loop, with a unique name.
        string fileout = filename;
        string argument = to_string(i); fileout.append("_10^"); fileout.append(argument); fileout.append(".txt");

        //Calculates the value of the steplength h.
        //And the value of h^2, to make the code more readable later on.
        double h = 1.0 / n;
        double hh = h * h;
       
        //Declaring arrays needed for calculations. Size of n+1 to accomodate known endpoints.
        //d is the vector for the diagonal of matrix A. x will contain all the values of x. b will contain all f_i(x_i).
        //solution will be filled with the calculated solution for different values of x.  
        double* d_array = new double[n + 1];
        double* x = new double[n + 1];
        double* solution = new double[n + 1];
        double* b_array = new double[n + 1];

        //Filling d and solutiion with known endpoints.
        d_array[0] = d_array[n] = 2; solution[0] = solution[n] = 0.0;

        //Fetches current time. Used in determining  execution time.
        auto start = chrono::high_resolution_clock::now();

        //Updates diagonal values (d) and calculates startingvalues of b.
        for (int i = 1; i < n; i++) {
            d_array[i] = (i + 1.0) / ((double)i);
        }
        for (int i = 0; i <= n; i++) {
            x[i] = i * h;
            b_array[i] = hh * f(i * h);
        }

        //Forward substitution
        //Updates all the values of the b-vector.
        for (int i = 2; i < n; i++) {
            b_array[i] = b_array[i] +  b_array[i - 1]/ d_array[i-1];
        }

        //Backward substitution
        //Calculates the solution for all x_i, and stores in array solution.
        solution[n - 1] = b_array[n - 1] / d_array[n - 1];
        for (int i = n - 2; i > 0; i--) solution[i] = (b_array[i] + solution[i+1])/d_array[i];
        
        //Fetches current time after algorithm is done. Then calculates execution time in microseconds.
        auto finish = chrono::high_resolution_clock::now();
        auto timeused = chrono::duration_cast<chrono::microseconds>(finish - start).count();

        //Updates execution time of  loop of power to 10 of 1.
        ExecutionTimes[i - 1] = timeused;

        //Variable to determine max relative error. Will always be larger than 0, so chooses initial value to be 0.
        double max_RE = 0.0;
        
        //Array for storing relativ error for all solutions.
        double *RelativeError = new double[n];

        //Calculates the relative error for all solutions. Stores them in the array RelativeError.
        //Updates the value of max_RE if it is larger than curren value.
        for (int i = 1; i < n; i++) {
            double xval = x[i];
            double RE = fabs((f_exact(xval) - solution[i]) / f_exact(xval));
            RelativeError[i - 1] = RE;

            if (RE > max_RE) max_RE = RE;
        }

        //Updates MaxRelativeError loop of power of 10 of i.
        MaxRelativeError[i - 1] = max_RE;

        //Opens file with desired name, and fills inn header.
        ofile.open(fileout);
        ofile << "x:                 approx:        exact:          relative error (log10)" << endl;
        
        //Writes all x_i, with corresponding calculated solution, exact solution and log10 of relative error.
        for (int i = 1; i < n; i++) {
            double xval = x[i];
            ofile << setw(15) << setprecision(8) << x[i];
            ofile << setw(15) << setprecision(8) << solution[i];
            ofile << setw(15) << setprecision(8) << f_exact(xval);
            ofile << setw(15) << setprecision(8) << log10(RelativeError[i-1]) << endl;
        }

        //Cloeses file.
        ofile.close();

        //Deletes the arrays, as to free memory.
        delete[] b_array; delete[] d_array; delete[] x; delete[] solution;
    }

    //Makes a file with name filename_executiontimes_maxRE.txt, 
    //and writes execution time and max relative error for different 10^n to file.
    filename.append("_executiontimes_maxRE"); filename.append(".txt");
    ofile.open(filename);
    ofile << "Number of steps:         Execution time [microseconds]:      Max relative error (log_10): " << endl;
    for (int i = 1; i <= exponent; i++) {
        ofile << setw(15) << setprecision(8) << pow(10.0, i);
        ofile << setw(30) << setprecision(8) << ExecutionTimes[i - 1];
        ofile << setw(28) << setprecision(8) << log10(MaxRelativeError[i - 1]) << endl;
    }
    ofile.close();

    return 0;
}

//Functions to calculate f_i. 
//Used in updating array b__array.
double f(double x)
{
    return 100.0 * exp(-10.0 * x);
}

//Function to calculate the exact solution. 
double f_exact(double x)
{
    return 1.0 - (1.0 - exp(-10.0)) * x - exp(-10.0 * x);
}
