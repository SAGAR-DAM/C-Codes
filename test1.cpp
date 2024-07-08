#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

// Function to compute the value of a polynomial at a given point z
complex<double> evaluatePolynomial(const vector<double>& coefficients, complex<double> z) 
{
    complex<double> result = 0.0;
    int n = coefficients.size() - 1; // degree of polynomial

    for (int i = 0; i <= n; ++i) {
        result += coefficients[i] * pow(z, n - i);
    }

    return result;
}

// Function to apply the Durand-Kerner method to find all roots of a polynomial
vector<complex<double>> findRoots(const vector<double>& coefficients, double epsilon = 1e-8, int maxIterations = 1000) 
{
    int n = coefficients.size() - 1; // degree of polynomial
    vector<complex<double>> roots(n);

    // Initial guesses for roots (usually on the unit circle in the complex plane)
    for (int i = 0; i < n; ++i) {
        roots[i] = polar(1.0, 2.0 * M_PI * i / n);
    }

    // Iterative refinement using Durand-Kerner method
    for (int iter = 0; iter < maxIterations; ++iter) {
        vector<complex<double>> newRoots(n);

        for (int i = 0; i < n; ++i) {
            complex<double> numerator = evaluatePolynomial(coefficients, roots[i]);
            complex<double> denominator = 1.0;

            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    denominator *= (roots[i] - roots[j]);
                }
            }

            newRoots[i] = roots[i] - numerator / denominator;
        }

        // Check convergence
        bool converged = true;
        for (int i = 0; i < n; ++i) {
            if (abs(newRoots[i] - roots[i]) > epsilon) {
                converged = false;
                break;
            }
        }

        if (converged) {
            return newRoots; // Roots found
        }

        roots = newRoots;
    }

    cerr << "Warning: Maximum iterations reached without convergence." << endl;
    return roots; // Return current roots (may not be accurate)
}

int main() 
{
    // Example usage: Solve x^3 - 6x^2 + 11x - 6 = 0
    vector<double> coefficients = { 4,0,1 }; // coefficients in descending order of powers

    vector<complex<double>> roots = findRoots(coefficients);

    cout << "Roots of the polynomial:" << endl;
    for (auto root : roots) {
        cout << root << endl;
    }

    // cout<<abs(std::complex<double>(1,1));

    return 0;
}
