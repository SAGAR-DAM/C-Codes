#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
#include <set>

// Function to evaluate the polynomial at a given complex number
std::complex<double> polynomial(const std::vector<double>& coeffs, std::complex<double> x) {
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


// Newton-Raphson method for complex numbers
std::complex<double> newtonRaphson(const std::vector<double>& coeffs, std::complex<double> initialSeed, int iterations) {
    std::complex<double> x = initialSeed;

    for (int i = 0; i < iterations; ++i) {
        std::complex<double> f_x = polynomial(coeffs, x);
        std::complex<double> f_prime_x = derivative(coeffs, x);

        if (f_prime_x == 0.0) {
            std::cerr << "Derivative is zero. No solution found." << std::endl;
            return x;
        }

        // Update x using Newton-Raphson formula
        x = x - f_x / f_prime_x;

        // Display the current iteration and value
        std::cout << "Iteration " << i + 1 << ": x = " << std::setprecision(10) << x << std::endl;
    }


    return x;
}



int main() {
    // Define the polynomial coefficients for x^5 + 3x^4 - 5x^3 + 10x^2 + 17x - 30
    std::vector<double> coefficients = {1,3,-5,10,17,-30}; // a0, a1, a2, a3, a4, a5 (for x^5)

    // User input for complex initial seed
    double realPart, imagPart;
    std::cout << "Enter the real part of the initial seed: ";
    std::cin >> realPart;
    std::cout << "Enter the imaginary part of the initial seed: ";
    std::cin >> imagPart;
    std::complex<double> initialSeed(realPart, imagPart);
    
    int iterations = 100; // Number of iterations

    // Apply Newton-Raphson method
    std::complex<double> root = newtonRaphson(coefficients, initialSeed, iterations);
    
    // Calculate f(root) after iterations
    std::complex<double> f_root = polynomial(coefficients, root);
    
    std::cout << "Approximate root: " << std::setprecision(10) << root << std::endl;
    std::cout << "Value of f(root): " << std::setprecision(10) << f_root << std::endl;

    return 0;
}
