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

// for lapack installation: https://chatgpt.com/share/6700c6a9-1158-8013-90c8-7e4225bd2664

/*
If lapack is installed in the location: "D:\C++modules\lapack-3.12.0"
and the nonempty folder: "D:\C++modules\lapack-3.12.0\build\lib" exists
Run the code with:  g++ -o .\polynomial_solver "D:\C++ codes\My_C++_modules\np.cpp"  .\polynomial_solver.cpp -LD:\C++modules\lapack-3.12.0\build\lib -llapack -lblas -lgfortran 
*/

extern "C" 
{
    void dgeev_(char*, char*, int*, double*, int*, double*, double*, double*, int*, double*, int*, double*, int*, int*);
}

using namespace std;

// Function to create a companion matrix for the polynomial
std::vector<double> createCompanionMatrix(const std::vector<double>& coefficients) 
{
    int n = coefficients.size() - 1;  // Degree of the polynomial

    if (coefficients[n] == 0) 
    {
        throw std::runtime_error("Error: Leading coefficient cannot be zero.");
    }

    // Companion matrix will be of size n x n
    std::vector<double> companionMatrix(n * n, 0.0);

    // Fill the last column with -coefficients / a_n (a_n is the leading coefficient)
    for (int i = 0; i < n; ++i) 
    {
        companionMatrix[i * n + (n - 1)] = -coefficients[i] / coefficients[n];  // Here we divide by the leading coefficient
    }

    // Fill the sub-diagonal with 1s
    for (int i = 1; i < n; ++i) 
    {
        companionMatrix[i * n + (i - 1)] = 1.0;
    }

    return companionMatrix;
}



// Function to compute roots of a polynomial using LAPACK
std::vector<std::complex<double>> solvePolynomial(const std::vector<double>& coefficients) 
{
    int degree = coefficients.size() - 1;

    if (degree == 0) 
    {
        throw std::runtime_error("Error: Not a valid polynomial (degree 0).");
    }

    std::vector<double> companionMatrix = createCompanionMatrix(coefficients);

    // Arrays for eigenvalues
    std::vector<double> wr(degree), wi(degree);

    // LAPACK parameters
    int n = degree;
    int lda = n;
    int info;
    int lwork = 4 * n;
    std::vector<double> work(4 * n);  // Work array

    // Temporary arrays for left and right eigenvectors (not needed)
    std::vector<double> vl(1), vr(1);
    char jobvl = 'N';  // No left eigenvectors are computed
    char jobvr = 'N';  // No right eigenvectors are computed

    // Call LAPACK routine to compute the eigenvalues (roots)
    dgeev_(&jobvl, &jobvr, &n, companionMatrix.data(), &lda, wr.data(), wi.data(), vl.data(), &lda, vr.data(), &lda, work.data(), &lwork, &info);

    // Check for success
    if (info != 0) 
    {
        throw std::runtime_error("Error: dgeev failed with info = " + std::to_string(info));
    }

    // Collect the inverted roots as complex numbers
    std::vector<std::complex<double>> roots;
    for (int i = 0; i < degree; ++i) 
    {
        // Invert the root (1 / root), considering both real and imaginary parts
        std::complex<double> root(wr[i], wi[i]);

        // Ensure we are not dividing by zero
        if (root == std::complex<double>(0.0, 0.0)) 
        {
            throw std::runtime_error("Error: Encountered zero root, cannot invert.");
        }

        // Invert the root
        roots.push_back(1.0 / root);
    }

    return roots;
}



// Function to generate a random double vector of a given size and range
std::vector<double> generateRandomDoubleVector(int size, double min_value, double max_value) {
    // Seed with a value that changes, like the current time in nanoseconds
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed); // Mersenne Twister engine

    // Define a uniform distribution for double values in the range [min_value, max_value]
    std::uniform_real_distribution<> distrib(min_value, max_value);

    // Create a vector of random doubles
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


// Function to convert a vector of integers to a vector of doubles
std::vector<double> convertIntToDouble(const std::vector<int>& int_vector) {
    std::vector<double> double_vector(int_vector.size());

    // Iterate through the integer vector and cast each element to double
    for (size_t i = 0; i < int_vector.size(); ++i) {
        double_vector[i] = static_cast<double>(int_vector[i]);
    }

    return double_vector;
}


// Function to construct the polynomial string from coefficients
std::string constructPolynomial(const std::vector<double>& coefficients) {
    std::ostringstream polynomial; // To hold the polynomial string
    int n = coefficients.size(); // Degree of the polynomial

    for (int i = 0; i < n; ++i) {
        double coeff = coefficients[i];
        int power = n - 1 - i; // Calculate the power

        // Skip zero coefficients
        if (coeff == 0) {
            continue;
        }

        // Add the coefficient and variable part
        if (polynomial.str().length() > 0) { // If not the first term
            if (coeff > 0) {
                polynomial << "+";
            } else {
                polynomial << "-";
            }
        }

        // Handle the absolute value of the coefficient
        if (abs(coeff) != 1 || power == 0) {
            polynomial << abs(coeff); // Only print coefficient if it's not ±1 or if it's constant
        }

        // Handle the variable and power part
        if (power > 0) {
            polynomial << "x"; // Add variable x
            if (power > 1) {
                polynomial << "^" << power; // Add power if it's greater than 1
            }
        }
    }

    return polynomial.str(); // Return the constructed polynomial as a string
}



// Function to create a Gnuplot script to plot the roots
void createGnuplotScript(const std::string& polynomial, const std::vector<std::complex<double>>& roots) {
    double min_real = std::numeric_limits<double>::max();
    double max_real = std::numeric_limits<double>::min();
    double min_imag = std::numeric_limits<double>::max();
    double max_imag = std::numeric_limits<double>::min();

    // Calculate min and max for real and imaginary parts
    for (const auto& root : roots) {
        double real_part = root.real();
        double imag_part = root.imag();
        
        if (real_part < min_real) min_real = real_part;
        if (real_part > max_real) max_real = real_part;
        if (imag_part < min_imag) min_imag = imag_part;
        if (imag_part > max_imag) max_imag = imag_part;
    }

    // Prepare polynomial title with proper formatting
    std::string formattedPolynomial = polynomial;
    
    // Use regex to format the polynomial string
    std::regex re("(\\w)\\^([0-9]+)"); // Matches variable^exponent
    formattedPolynomial = std::regex_replace(formattedPolynomial, re, "$1^{ $2 }"); // Replace with correct formatting

    std::ostringstream script;
    script << "set title \"" << formattedPolynomial << " = 0"<<"\"\n";
    script << "set xlabel \"Real\"\n";
    script << "set ylabel \"Imaginary\"\n";
    script << "set grid\n";
    
    // Set dynamic ranges based on min and max values
    script << "set xrange [" << 1.1*min_real << ":" << 1.1*max_real << "]\n"; // Add some padding
    script << "set yrange [" << 1.1*min_imag << ":" << 1.1*max_imag << "]\n"; // Add some padding
    script << "plot '-' with points pt 7 ps 1.5 title 'Roots'\n";

    // Add the roots to the Gnuplot script
    for (const auto& root : roots) {
        script << root.real() << " " << root.imag() << "\n";
    }
    script << "e\n"; // End of data
    
    // Add pause to keep the plot window open
    script << "pause -1\n"; // Pause indefinitely until user closes the window

    // Write the script to a file
    std::ofstream gnuplotScript("plot_roots.gp");
    gnuplotScript << script.str();
    gnuplotScript.close();
}




int main() 
{
    // Example: Solve the polynomial x^4 - 1 = 0
    // std::vector<double> coefficients = {1,0,3,0,4}; // Represents x^4 - 1
    int size;
    cout<<"Enter size: ";
    cin>>size;
    cout<<endl;
    std::vector<double> coefficients = convertIntToDouble(generateRandomVector(size+1,-10,10));

    // Construct the polynomial string
    std::string polynomial = constructPolynomial(coefficients);

    // Print the polynomial equation
    std::cout << "The equation:\n";
    std::cout << polynomial << " = 0\n" << std::endl;
    try 
    {
        std::vector<std::complex<double>> roots = solvePolynomial(coefficients);

        // Output the roots
        std::cout << "The roots of the polynomial are:\n";
        for (size_t i = 0; i < roots.size(); i++) 
        {
            std::cout << "Root " << i + 1 << ": " << roots[i] << std::endl;
        }

        // Create Gnuplot script to plot the roots
        createGnuplotScript(polynomial, roots);

        // Execute the Gnuplot script
        system("gnuplot plot_roots.gp");
    } 
    catch (const std::exception& e) 
    {
        std::cerr << e.what() << std::endl;
    }

    cout<<"press any key to exit"<<endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore leftover input
    std::cin.get();
    return 0;
}
