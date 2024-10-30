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



// Function to create a Gnuplot script to plot the roots
void createGnuplotScript(const std::string& polynomial, const std::vector<std::complex<double>>& roots) 
{
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
    script << "set title \"" << "Eigenvalues of the given matrix (roots of characteristic eq)\\n" << formattedPolynomial << " = 0"<<"\"\n";
    script << "set xlabel \"Real\"\n";
    script << "set ylabel \"Imaginary\"\n";
    script << "set grid\n";
    
    // Set dynamic ranges based on min and max values
    script << "set xrange [" << -0.5+1.1*min_real << ":" << 0.5+1.1*max_real << "]\n"; // Add some padding
    script << "set yrange [" << -0.5+1.1*min_imag << ":" << 0.5+1.1*max_imag << "]\n"; // Add some padding
    script << "plot '-' with points pt 7 ps 1.5 lc rgb 'red' title 'Roots'\n";


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
    int size;
    cout<<"Enter size: ";
    cin>>size;
    cout<<endl;
    
    Matrix A = getRandomIntMatrix(size,size,-10,10);
    std::cout << "\033[31m"; // Set text color to red
    cout<<"========================================================================================================================"<<endl;
    cout<<"The Matrix:"<<endl<<"-------------"<<endl<<A<<endl<<endl;
    cout<<"========================================================================================================================"<<endl<<endl;

    //get the characteristic polynomial
    std::vector<double> coefficients = Matrix::coeff_charac(A);

    // Construct the polynomial string
    std::string polynomial = constructPolynomial(coefficients);

    // Print the polynomial equation
    std::cout << "The Characteristic Polynomial:\n";
    std::cout << polynomial << " = 0\n" << std::endl;

    try 
    {
        int seed_number = 125;
        if(seed_number>5)
        {
            seed_number=size*size*size;
        }
        
        vector<complex<double>> seeds = create_circ_root(seed_number);
        std::vector<std::complex<double>> roots; // Vector to store the final roots
        int max_iteration=1000;
        if(size>10)
        {
            max_iteration = size*size*size;
        }

        int counter = 0;
        // Loop over each seed to find a root using predictRoot
        for (const auto& seed : seeds) 
        {
            std::complex<double> root = predictRoot(coefficients, seed, 1e-10, max_iteration); // tolerance = 1e-6, maxIterations = 100
            if(abs(evalpolynomial(coefficients,root))<0.1)
            {
                roots.push_back(root);
            }
            counter = counter+1 ;
            if(removeCloseDuplicates(roots).size()==size)
            {   
                cout<<"Loop run counter: "<<counter<<endl<<endl;
                std::cout << "\033[32m"; // Set text color to green
                cout<<"========================================================================================================================"<<endl<<endl;
                break;
            }
        }

        // Remove duplicates within a certain tolerance
        std::vector<std::complex<double>> uniqueRoots = removeCloseDuplicates(roots);

        for (size_t i = 0; i < uniqueRoots.size(); ++i) {
            std::complex<double> root = uniqueRoots[i];
            
            // Calculate the error if required (assuming a function characPolyError exists)
            double error = std::abs(evalpolynomial(coefficients, root));  // Example error calculation

            // Print the root index, root value, and error, all justified
            std::cout << "\033[34m" << std::left << "Root " << i + 1 << ": " << "\033[32m"
                    << std::setw(1) << std::fixed << std::setprecision(8) << "(" << root.real() << ", " << root.imag() << ")"
                    << "\033[33m" << "\twith Error Charac_Poly(x) ~ " << std::scientific << std::setprecision(2) << error << "\n";
        }
        
        cout<<"\033[32m"<<"========================================================================================================================"<<endl;
        // Reset the color back to default
        std::cout << "\033[0m";

        // Create Gnuplot script to plot the unique roots
        createGnuplotScript(polynomial, uniqueRoots);

        // Execute the Gnuplot script
        system("gnuplot plot_roots.gp");
    } 
    catch (const std::exception& e) 
    {
        cout<<"error occured !!"<<endl;
        std::cerr << e.what() << std::endl;
    }


    return 0;
}