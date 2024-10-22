#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <regex>

using namespace std;

// for lapack installation: https://chatgpt.com/share/6700c6a9-1158-8013-90c8-7e4225bd2664

/*
If lapack is installed in the location: "D:\C++modules\lapack-3.12.0"
and the nonempty folder: "D:\C++modules\lapack-3.12.0\build\lib" exists
Run the code with:  g++ -o .\polynomial_solver "D:\C++ codes\My_C++_modules\np.cpp"  .\polynomial_solver.cpp -LD:\C++modules\lapack-3.12.0\build\lib -llapack -lblas -lgfortran 
*/


extern "C" {
    // LAPACK routine to solve a system of linear equations Ax = B
    void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
}

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



// Function to generate a random system of linear equations
std::pair<std::vector<std::vector<double>>, std::vector<double>> generateRandomEquations(int num_variables) {
    std::vector<std::vector<double>> A(num_variables, std::vector<double>(num_variables));
    std::vector<double> B = convertIntToDoubleVector(generateRandomVector(num_variables, -10.0*num_variables, 15.0*num_variables));

    A = convertIntMatrixToDouble(generateRandomIntMatrix(num_variables,num_variables,-20,20));

    return {A, B};
}

// Function to print the system of equations
void printSystem(const std::vector<std::vector<double>>& A, const std::vector<double>& B) {
    int n = A.size();
    std::ostringstream equation;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (A[j][i] >= 0 && j > 0) {
                equation << "+";
            }
            equation << A[j][i] << "x" << (j + 1) << " ";
        }
        equation << "= " << B[i] << std::endl;
    }

    std::cout << "System of equations:\n";
    std::cout << equation.str() << std::endl;
}

// Function to solve the system of linear equations using LAPACK
std::vector<double> solveLinearSystem(std::vector<std::vector<double>>& A, std::vector<double>& B) {
    int n = A.size(); // Number of variables
    int nrhs = 1; // Number of right-hand sides (we have only one)
    int lda = n;
    int ldb = n;
    int info;
    std::vector<int> ipiv(n); // Pivot indices

    // Convert the matrix A to a 1D array
    std::vector<double> A_flat(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A_flat[i * n + j] = A[i][j];
        }
    }

    // Call LAPACK to solve the system
    dgesv_(&n, &nrhs, A_flat.data(), &lda, ipiv.data(), B.data(), &ldb, &info);

    if (info != 0) {
        throw std::runtime_error("Error: dgesv failed with info = " + std::to_string(info));
    }

    return B; // B now contains the solution
}

int main() {
    int num_variables;
    std::cout << "Enter the number of variables: ";
    std::cin >> num_variables;

    std::cout << "Do you want to input the system of equations manually? (y/n): ";
    char choice;
    std::cin >> choice;

    std::vector<std::vector<double>> A;
    std::vector<double> B;

    if (choice == 'y' || choice == 'Y') {
        A.resize(num_variables, std::vector<double>(num_variables));
        B.resize(num_variables);
        std::cout << "Enter the coefficients of the system of equations:" << std::endl;
        for (int i = 0; i < num_variables; ++i) {
            std::cout << "Equation " << (i + 1) << " (Enter coefficients for x1, x2,... x" << num_variables << " and constant term):" << std::endl;
            for (int j = 0; j < num_variables; ++j) {
                std::cin >> A[i][j];
            }
            std::cin >> B[i];
        }
    } else {
        std::tie(A, B) = generateRandomEquations(num_variables);
    }

    printSystem(A, B);

    try {
        std::vector<double> solution = solveLinearSystem(A, B);

        std::cout << "The solution to the system is:" << std::endl;
        for (int i = 0; i < solution.size(); ++i) {
            std::cout << "x" << (i + 1) << " = " << solution[i] << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }

    std::cout<<"Press any key to exit";

    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore leftover input
    std::cin.get();
    return 0;
}
