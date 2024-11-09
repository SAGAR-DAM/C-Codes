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
#include <stdexcept>
#include <regex>
#include <iomanip>
#include <chrono>   // to create new random numbers each time. Visit: https://chatgpt.com/share/6700cacc-8640-8013-bc41-94054ab0dc0b

using namespace std;


const int ITERATION = 50;
const int RES = 700;


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



// Function to create n-th roots of unity
vector<complex<double>> create_circ_root(int n) 
{
    vector<complex<double>> roots;
    double angle_increment = 2.0 * M_PI / n;

    for (int k = 0; k < n; ++k) 
    {
        double angle = k * angle_increment;
        complex<double> z(50*cos(angle), 50*sin(angle));
        roots.push_back(z);
    }

    return roots;
}



// Function to construct the polynomial string from coefficients
std::string constructPolynomial(const std::vector<double>& coefficients) 
{
    std::ostringstream polynomial;
    int n = coefficients.size();

    for (int i = 0; i < n; ++i) {
        double coeff = coefficients[i];
        int power = n - 1 - i;

        // Skip zero coefficients
        if (coeff == 0) {
            continue;
        }

        // Handle the sign for the first term or subsequent terms
        if (polynomial.str().empty()) { // First term
            if (coeff < 0) {
                polynomial << "-";
            }
        } else { // Subsequent terms
            polynomial << (coeff > 0 ? "+" : "-");
        }

        // Handle the absolute value of the coefficient
        if (abs(coeff) != 1 || power == 0) {
            polynomial << abs(coeff);
        }

        // Handle the variable and power part
        if (power > 0) {
            polynomial << "x";
            if (power > 1) {
                polynomial << "^" << power;
            }
        }
    }
    return polynomial.str();
}



// Function to evaluate the polynomial at a given complex number
std::complex<double> evalpolynomial(const std::vector<double>& coeffs, std::complex<double> x) {
    std::complex<double> result = 0.0;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * std::pow(x, coeffs.size() - 1 - i);
    }
    return result;
}



// Function to evaluate the derivative of the polynomial at a given complex number
std::complex<double> derivative(const std::vector<double>& coeffs, std::complex<double> x) {
    std::complex<double> result = 0.0;
    for (size_t i = 0; i < coeffs.size() - 1; ++i) {
        result += coeffs[i] * static_cast<std::complex<double>>(coeffs.size() - 1 - i) * std::pow(x, coeffs.size() - 2 - i);
    }
    return result;
}



// Newton-Raphson method with tolerance to predict root
std::complex<double> predictRoot(const std::vector<double>& coeffs, std::complex<double> initialSeed, double tolerance, int maxIterations) {
    std::complex<double> x = initialSeed;

    for (int i = 0; i < maxIterations; ++i) {
        std::complex<double> f_x = evalpolynomial(coeffs, x);
        std::complex<double> f_prime_x = derivative(coeffs, x);

        if (f_prime_x == 0.0) {
            std::cerr << "Derivative is zero. No solution found." << std::endl;
            return x;
        }

        // Update x using Newton-Raphson formula
        x = x - f_x / f_prime_x;

        // Check if the error (|f(x)|) is within tolerance
        if (std::abs(f_x) < tolerance) {
            //std::cout << "Converged in " << i + 1 << " iterations." << std::endl;
            break;
        }

        // Display the current iteration and value
        //std::cout << "Iteration " << i + 1 << ": x = " << std::setprecision(10) << x << ", |f(x)| = " << std::abs(f_x) << std::endl;
    }

    return x;
}



// Function to remove duplicate roots that are within a specified tolerance
std::vector<std::complex<double>> removeCloseDuplicates(const std::vector<std::complex<double>>& roots, double tolerance = 1e-10) 
{
    std::vector<std::complex<double>> uniqueRoots;

    for (const auto& root : roots) {
        bool isDuplicate = false;
        
        // Check against all unique roots found so far
        for (auto& uniqueRoot : uniqueRoots) {
            if (std::abs(root - uniqueRoot) < tolerance) {
                // If root is close to uniqueRoot, average them and mark as duplicate
                uniqueRoot = (uniqueRoot + root) / 2.0;
                isDuplicate = true;
                break;
            }
        }

        // If no close duplicate found, add root to uniqueRoots
        if (!isDuplicate) {
            uniqueRoots.push_back(root);
        }
    }

    // Set small imaginary parts to zero for all unique roots
    for (auto& root : uniqueRoots) {
        if (std::abs(root.imag()) < tolerance) {
            root = std::complex<double>(root.real(), 0.0);  // Set imaginary part to 0
        }
    }

    return uniqueRoots;
}



// Function to compute f(z)
complex<double> f(complex<double>& z, const vector<complex<double>>& root) {
    complex<double> val = 1.0;
    for (const auto& r : root) {
        val *= (z - r);
    }
    return val;
}



// Function to compute df(z)
complex<double> df(complex<double>& z, const vector<complex<double>>& root) {
    complex<double> val = 0.0;
    for (size_t i = 0; i < root.size(); ++i) 
    {
        complex<double> mult = 1.0;
        for (size_t j = 0; j < root.size(); ++j) 
        {
            if (j != i) {
                mult *= (z - root[j]);
            }
        }
        val += mult;
    }
    return val;
}



// Function to perform a single iteration step
complex<double> iteration_step(complex<double>& z, const vector<complex<double>>& root) 
{
    if(abs(f(z,root))<=0.000001 || abs(df(z,root))<=0.000001)
    {
        return z;
    }
    return z - f(z, root) / df(z, root);
}




int main() 
{
    int size;
    cout << "\033[34m" << "Enter order of equation: ";
    cin >> size;

    cout << endl;
    vector<double> coefficients = convertIntToDouble(generateRandomVector(size + 1, -10, 10));
    if (coefficients[0] == 0) 
    {
        coefficients[0] = 1; // Prevent zero dominant coefficient
    }

    // Construct and display the polynomial equation
    string polynomial = constructPolynomial(coefficients);
    cout << "\033[31m" << "The equation:\n" << polynomial << " = 0\n" << endl;

    vector<complex<double>> seeds = create_circ_root(size * size * size);
    vector<complex<double>> roots; // Store final roots

    int max_iteration = (size > 10) ? size * size * size : 1000;
    for (const auto& seed : seeds) 
    {
        complex<double> root = predictRoot(coefficients, seed, 1e-10, max_iteration);
        if (abs(evalpolynomial(coefficients, root)) < 0.1) 
        {
            roots.push_back(root);
        }
        if (removeCloseDuplicates(roots).size() == size) 
        {
            cout << "Roots found within loop limit.\n";
            break;
        }
    }

    vector<complex<double>> uniqueRoots = removeCloseDuplicates(roots);
    cout << "\033[34mRoots of the polynomial equation:\n---------------------------------\n\n"<<"\033[32m";
    for (size_t i = 0; i < uniqueRoots.size(); ++i) 
    {
        complex<double> root = uniqueRoots[i];
        double error = abs(evalpolynomial(coefficients, root));
        cout << "Root " << i + 1 << ": " << fixed << setprecision(8) 
             << "(" << root.real() << ", " << root.imag() << ")"
             << "\tError ~ " << scientific << setprecision(2) << error << "\n";
    }
    cout << "\033[0m";

    // Plot Newton Fractal using root data
    int rootnumber = uniqueRoots.size();
    double maxReal = 1.1 * max(abs(real(uniqueRoots[0])), abs(imag(uniqueRoots[0])));
    for (const auto& r : uniqueRoots) 
    {
        maxReal = max(maxReal, 1.1 * max(abs(real(r)), abs(imag(r))));
    }


    vector<double> x(2 * RES * (maxReal + 1));
    vector<double> y(2 * RES * (maxReal + 1));
    for (int i = 0; i < x.size(); ++i) 
    {
        x[i] = -maxReal + (2 * maxReal / (x.size() - 1)) * i;
    }
    for (int i = 0; i < y.size(); ++i) 
    {
        y[i] = -maxReal + (2 * maxReal / (y.size() - 1)) * i;
    }

    vector<complex<double>> z(x.size() * y.size());
    for (int i = 0; i < x.size(); i++) 
    {
        for (int j = 0; j < y.size(); j++) 
        {
            z[i * x.size() + j] = complex<double>(x[i], y[j]);
        }
    }


    for (int i = 0; i < ITERATION; i++) 
    {
        for (size_t j = 0; j < z.size(); j++) 
        {
            z[j] = iteration_step(z[j], uniqueRoots);
        }
    }

    vector<int> f_values(z.size());
    for (size_t i = 0; i < z.size(); i++) 
    {
        double dist = abs(z[i] - uniqueRoots[0]);
        int minindex = 0;
        for (int j = 1; j < uniqueRoots.size(); j++) 
        {
            double val = abs(z[i] - uniqueRoots[j]);
            if (val <= dist) 
            {
                minindex = j;
                dist = val;
            }
        }
        f_values[i] = minindex;
    }

    ofstream dataFile("data.txt");
    for (int i = 0; i < x.size(); i++) 
    {
        for (int j = 0; j < y.size(); j++) 
        {
            dataFile << x[i] << " " << y[j] << " " << f_values[i * x.size() + j] << "\n";
        }
    }
    dataFile.close();

    ofstream dataFile1("roots.txt");
    for (const auto& root : uniqueRoots) 
    {
        dataFile1 << real(root) << " " << imag(root) << " " << 1 << "\n";
    }
    dataFile1.close();

    ofstream gnuplotScript("plot.gnuplot");
    gnuplotScript << "set term wxt size 900,900\n"
                  << "set view map\n"
                  << "set size ratio -1\n"
                  << "set palette defined (0 0 0 0, 0.1 0 0 0.5, 0.2 0 0 1, 0.3 0 0.5 1, 0.4 0 1 1, "
                  << "0.5 0.5 1 0.5, 0.6 1 1 0, 0.7 1 0.5 0, 0.8 1 0 0, 0.9 0.5 0 0, 1 0.5 0 0)\n"
                  << "set xlabel 'Re(z)'\n"
                  << "set ylabel 'Im(z)'\n"
                  << "set xrange [" << -maxReal << ":" << maxReal << "]\n"
                  << "set yrange [" << -maxReal << ":" << maxReal << "]\n"
                  << "set title 'Newton Fractal'\n"
                  << "unset key\n"
                  << "splot 'data.txt' using 1:2:3 with image, 'roots.txt' using 1:2:(0) with points pt 7 ps 2 lc rgb 'white'\n";
    gnuplotScript.close();

    system("gnuplot -persist plot.gnuplot");

    cout << "\033[0m" << "Plot generated. Press enter to exit...\n";
    cin.get();

        // Delete the data files
    if (remove("data.txt") != 0) 
    {
        std::cerr << "Error deleting data.txt" << std::endl;
    } 
    else 
    {
        std::cout << "Deleted data.txt" << std::endl;
    }

    // Delete the data files
    if (remove("roots.txt") != 0) 
    {
        std::cerr << "Error deleting roots.txt" << std::endl;
    } 
    else 
    {
        std::cout << "Deleted roots.txt" << std::endl;
    }

    return 0;
}
