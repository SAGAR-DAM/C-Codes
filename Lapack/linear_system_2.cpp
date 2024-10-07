#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include "D:\C++ codes\My_C++_modules\np.hpp"
#include <limits>  
#include <regex>
#include <chrono>   // to create new random numbers each time. Visit: https://chatgpt.com/share/6700cacc-8640-8013-bc41-94054ab0dc0b

using namespace std;

// for lapack installation: https://chatgpt.com/share/6700c6a9-1158-8013-90c8-7e4225bd2664

/*
If lapack is installed in the location: "D:\C++modules\lapack-3.12.0"
and the nonempty folder: "D:\C++modules\lapack-3.12.0\build\lib" exists
Run the code with:  g++ -o .\polynomial_solver "D:\C++ codes\My_C++_modules\np.cpp"  .\polynomial_solver.cpp -LD:\C++modules\lapack-3.12.0\build\lib -llapack -lblas -lgfortran 
*/


// Function to generate a single random double in a given range
double generateRandomDouble(double min_value, double max_value) 
{
    // Seed with a value that changes, like the current time in nanoseconds
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed); // Mersenne Twister engine

    // Define a uniform distribution for double values in the range [min_value, max_value]
    std::uniform_real_distribution<> distrib(min_value, max_value);

    // Generate and return a single random double
    return distrib(gen);
}


// Function to generate a single random integer in a given range
int generateRandomNumber(int min_value, int max_value) 
{
    // Seed with a value that changes, like the current time in nanoseconds
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed); // Mersenne Twister engine

    // Define a distribution in the range [min_value, max_value]
    std::uniform_int_distribution<> distrib(min_value, max_value);

    // Generate and return a single random number
    return distrib(gen);
}


// Function to generate a random double vector
std::vector<double> generateRandomDoubleVector(int size, double min_value, double max_value) {
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed); // Mersenne Twister engine
    std::uniform_real_distribution<> distrib(min_value, max_value);

    std::vector<double> random_vector(size);
    for (int i = 0; i < size; ++i) {
        random_vector[i] = distrib(gen);
    }

    return random_vector;
}


// Function to generate a random integer vector of a given size and range
std::vector<int> generateRandomVector(int size, int min_value, int max_value) 
{
    // Seed with a value that changes, like the current time in nanoseconds
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed); // Mersenne Twister engine

    // Define a distribution in the range [min_value, max_value]
    std::uniform_int_distribution<> distrib(min_value, max_value);

    // Create a vector of random integers
    std::vector<int> random_vector(size);
    for (int i = 0; i < size; ++i) {
        random_vector[i] = distrib(gen);
    }

    return random_vector;
}



// Function to generate a random integer matrix (2D vector)
std::vector<std::vector<int>> generateRandomIntMatrix(int rows, int cols, int min_value, int max_value) 
{
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed); // Mersenne Twister engine
    std::uniform_int_distribution<> distrib(min_value, max_value);

    // Create a 2D vector (matrix) with the given number of rows and columns
    std::vector<std::vector<int>> random_matrix(rows, std::vector<int>(cols));

    // Fill the matrix with random integer values
    for (int i = 0; i < rows; ++i) 
    {
        for (int j = 0; j < cols; ++j) 
        {
            random_matrix[i][j] = distrib(gen); // Generate a random value for each element
        }
    }

    return random_matrix;
}



// Function to generate a random double matrix (2D vector)
std::vector<std::vector<double>> generateRandomDoubleMatrix(int rows, int cols, double min_value, double max_value) 
{
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed); // Mersenne Twister engine
    std::uniform_real_distribution<> distrib(min_value, max_value);

    // Create a 2D vector (matrix) with the given number of rows and columns
    std::vector<std::vector<double>> random_matrix(rows, std::vector<double>(cols));

    // Fill the matrix with random double values
    for (int i = 0; i < rows; ++i) 
    {
        for (int j = 0; j < cols; ++j) 
        {
            random_matrix[i][j] = distrib(gen); // Generate a random value for each element
        }
    }

    return random_matrix;
}



// Function to convert an integer to a double
double convertIntToDouble(int value) {
    // Use static_cast to convert the integer to a double
    return static_cast<double>(value);
}


// Function to convert a vector of integers to a vector of doubles
std::vector<double> convertIntToDoubleVector(const std::vector<int>& int_vector) {
    std::vector<double> double_vector(int_vector.size());

    // Iterate through the integer vector and cast each element to double
    for (size_t i = 0; i < int_vector.size(); ++i) {
        double_vector[i] = static_cast<double>(int_vector[i]);
    }

    return double_vector;
}


// Function to convert an integer matrix to a double matrix
std::vector<std::vector<double>> convertIntMatrixToDouble(const std::vector<std::vector<int>>& int_matrix) 
{
    int rows = int_matrix.size(); // Number of rows in the matrix
    int cols = int_matrix[0].size(); // Number of columns in the matrix

    // Create a 2D vector (matrix) for double values with the same size as the integer matrix
    std::vector<std::vector<double>> double_matrix(rows, std::vector<double>(cols));

    // Convert each element from int to double
    for (int i = 0; i < rows; ++i) 
    {
        for (int j = 0; j < cols; ++j) 
        {
            double_matrix[i][j] = static_cast<double>(int_matrix[i][j]);
        }
    }

    return double_matrix;
}



// Helper function to convert std::vector<std::vector<int>> to Matrix object
Matrix convertIntMatrixToMatrix(const std::vector<std::vector<int>>& int_matrix) {
    std::vector<std::vector<double>> double_matrix(int_matrix.size(), std::vector<double>(int_matrix[0].size()));

    // Convert each element from int to double
    for (int i = 0; i < int_matrix.size(); ++i) {
        for (int j = 0; j < int_matrix[0].size(); ++j) {
            double_matrix[i][j] = static_cast<double>(int_matrix[i][j]);
        }
    }

    // Create a Matrix object using the double matrix
    return Matrix(double_matrix);
}

// Function to generate a random system of linear equations
std::pair<Matrix, Matrix> generateRandomEquations(int num_variables) {
    // Generate random integer matrix for A (num_variables x num_variables)
    std::vector<std::vector<int>> random_int_matrix = generateRandomIntMatrix(num_variables, num_variables, -20, 20);

    // Convert integer matrix to double matrix and construct the Matrix object A
    Matrix A = convertIntMatrixToMatrix(random_int_matrix);

    // Generate random vector B (num_variables x 1) as a matrix
    std::vector<std::vector<int>> random_int_vector_B = generateRandomIntMatrix(num_variables, 1, -10.0 * num_variables, 15.0 * num_variables);

    // Convert the random integer vector B to Matrix
    Matrix B = convertIntMatrixToMatrix(random_int_vector_B);

    // Return both A and B as Matrix objects
    return {A, B};
}

// Function to print the system of equations
void printSystem(Matrix A, Matrix B) 
{
    int n = A.getRows();
    //std::ostringstream equation;
    std::cout << "System of equations:\n\n";

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (A[i][j] >= 0 && j > 0) {
                cout << "+";
            }
            cout << A[i][j] << "x" << (j + 1) << " ";
        }
        cout << " = " << B[i][0] << endl;
    }
    cout<<endl<<endl;
}



// Function to set elements less than epsilon to zero
void setSmallValuesToZero(Matrix& mat, double epsilon) {
    for (int i = 0; i < mat.getRows(); ++i) {
        for (int j = 0; j < mat.getCols(); ++j) {
            if (std::abs(mat[i][j]) < epsilon) {
                mat[i][j] = 0; // Set small values to zero
            }
        }
    }
}

int main()
{
    int num_variables;
    std::cout << "Enter the number of variables: ";
    std::cin >> num_variables;

    // Call the function to generate random equations and assign A and B separately
    std::pair<Matrix, Matrix> randomSystem = generateRandomEquations(num_variables);
    Matrix A = randomSystem.first;  // Coefficient matrix
    Matrix B = randomSystem.second; // Constants matrix (as a column vector)

    printSystem(A,B);

    Matrix inv_A = A.inverse();

    Matrix sol = Matrix::multiply(inv_A,B);

    for(int i=0;i<num_variables;++i)
    {
        cout<<"x"<<(i+1)<<" = "<<sol[i][0]<<endl;
    }

    Matrix cross_varify = Matrix::subtract(Matrix::multiply(A,sol),B).transpose();
    double error_range=0;
    for(int i=0; i<num_variables; ++i)
    {
        error_range = error_range + abs(cross_varify[0][i]);
    }
    setSmallValuesToZero(cross_varify, 1e-8); // Call the function to set small values to zero

    cout<<endl<<"cross_varify: A.sol = "<<"("<<cross_varify<<")"<<endl;
    cout<<"With error level: "<<error_range/num_variables;
    return 0;
}
