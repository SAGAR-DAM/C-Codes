#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <chrono>
#include <stdexcept>
#include <cmath>
#include <iomanip>

// se: https://pages.hmc.edu/ruye/MachineLearning/lectures/ch2/node7.html

using namespace std;

// Structure to represent a complex vector
using ComplexVector = std::vector<std::complex<double>>;

// Function to perform Gaussian elimination
ComplexVector gaussianElimination(const std::vector<std::vector<std::complex<double>>>& A, const ComplexVector& b) {
    int n = A.size();
    std::vector<std::vector<std::complex<double>>> augmented(n, std::vector<std::complex<double>>(n + 1));

    // Fill the augmented matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][n] = b[i]; // Last column is vector b
    }

    // Forward elimination
    for (int i = 0; i < n; ++i) {
        // Find the pivot row
        for (int k = i + 1; k < n; ++k) {
            if (abs(augmented[k][i]) > abs(augmented[i][i])) {
                swap(augmented[i], augmented[k]);
            }
        }

        // Make the elements below the pivot 0
        for (int k = i + 1; k < n; ++k) {
            std::complex<double> factor = augmented[k][i] / augmented[i][i];
            for (int j = i; j < n + 1; ++j) {
                augmented[k][j] -= factor * augmented[i][j];
            }
        }
    }

    // Back substitution
    ComplexVector x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = augmented[i][n] / augmented[i][i];
        for (int j = i - 1; j >= 0; --j) {
            augmented[j][n] -= augmented[j][i] * x[i];
        }
    }

    return x;
}

// Function to generate a vector of complex numbers with real and imaginary parts in range [-20, 20]
std::vector<std::complex<double>> generateComplexInitialSeeds(int num_vars, int min_value = -20, int max_value = 20) {
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> distrib(min_value, max_value);
    std::vector<std::complex<double>> initial_seeds(num_vars);

    for (int i = 0; i < num_vars; ++i) {
        initial_seeds[i] = {distrib(gen), distrib(gen)}; // Real and Imaginary parts both in range [-20, 20]
    }
    return initial_seeds;
}

// Function to generate a random polynomial system for multiple variables with integer coefficients
std::vector<std::vector<std::vector<double>>> generateRandomPolynomialSystem(int num_vars, int max_degree) {
    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<> distrib(-10, 10); // Range for integer coefficients

    std::vector<std::vector<std::vector<double>>> system(num_vars);
    for (int i = 0; i < num_vars; ++i) {
        system[i].resize(num_vars);
        for (int j = 0; j < num_vars; ++j) {
            system[i][j].resize(max_degree + 1);
            for (int k = 0; k <= max_degree; ++k) {
                system[i][j][k] = static_cast<double>(distrib(gen)); // Store as double
            }
        }
    }
    return system;
}

// Function to evaluate a polynomial equation for multiple complex variables
std::complex<double> evaluatePolynomial(const std::vector<std::vector<double>>& equation, const std::vector<std::complex<double>>& variables) {
    std::complex<double> result = 0.0;
    for (int var_index = 0; var_index < equation.size(); ++var_index) {
        int degree = equation[var_index].size() - 1;
        for (int d = 0; d <= degree; ++d) {
            result += equation[var_index][d] * pow(variables[var_index], degree - d);
        }
    }
    return result;
}

// Compute partial derivatives for each variable at a given point for Newton's method
std::vector<std::vector<std::complex<double>>> computeJacobian(const std::vector<std::vector<std::vector<double>>>& system, const std::vector<std::complex<double>>& variables) {
    int num_vars = system.size();
    std::vector<std::vector<std::complex<double>>> jacobian(num_vars, std::vector<std::complex<double>>(num_vars));

    for (int i = 0; i < num_vars; ++i) {
        for (int j = 0; j < num_vars; ++j) {
            std::complex<double> h = 1e-8;
            std::vector<std::complex<double>> variables_forward = variables;
            variables_forward[j] += h;
            std::complex<double> f_forward = evaluatePolynomial(system[i], variables_forward);
            std::complex<double> f_current = evaluatePolynomial(system[i], variables);
            jacobian[i][j] = (f_forward - f_current) / h;
        }
    }
    return jacobian;
}

// Newton-Raphson Method for solving polynomial system in complex numbers
std::vector<std::complex<double>> newtonsMethod(
    const std::vector<std::vector<std::vector<double>>>& system,
    std::vector<std::complex<double>> initial_guess,
    double tolerance,
    int maxIterations
) {
    int num_vars = system.size();
    std::vector<std::complex<double>> variables = initial_guess;

    for (int iter = 0; iter < maxIterations; ++iter) {
        // Evaluate each polynomial at the current guess
        std::vector<std::complex<double>> f_values(num_vars);
        for (int i = 0; i < num_vars; ++i) {
            f_values[i] = evaluatePolynomial(system[i], variables);
        }

        // Check if the solution meets the tolerance
        double max_error = 0.0;
        for (const auto& val : f_values) {
            max_error = max(max_error, abs(val));
        }
        if (max_error < tolerance) {
            cout << "Converged in " << iter + 1 << " iterations." << endl;
            return variables;
        }

        // Compute Jacobian matrix at the current guess
        std::vector<std::vector<std::complex<double>>> jacobian = computeJacobian(system, variables);

        // Solve the linear system using Gaussian elimination
        std::vector<std::complex<double>> delta = gaussianElimination(jacobian, f_values);

        // Update guess
        for (int i = 0; i < num_vars; ++i) {
            variables[i] += delta[i];
        }
    }
    throw std::runtime_error("Newton's method did not converge");
}

int main() {
    int num_vars, max_degree;
    cout << "Enter number of variables: ";
    cin >> num_vars;
    cout << "Enter maximum degree of the polynomial: ";
    cin >> max_degree;

    // Generate a random polynomial system with integer coefficients
    std::vector<std::vector<std::vector<double>>> system = generateRandomPolynomialSystem(num_vars, max_degree);

    cout << "\nGenerated Polynomial System:\n";
    for (int i = 0; i < num_vars; ++i) {
        cout << "Equation " << i + 1 << ": ";
        double constant_term_sum = 0.0;  // To accumulate all constant terms
        
        for (int j = 0; j < num_vars; ++j) {
            char variable_name = 'a' + j;
            for (int k = 0; k < system[i][j].size(); ++k) {
                int degree = max_degree - k;
                if (degree == 0) {
                    constant_term_sum += system[i][j][k];
                } else {
                    if (system[i][j][k] >= 0 && (j > 0 || k > 0)) cout << "+";
                    cout << system[i][j][k] << variable_name;
                    if (degree > 1) cout << "^" << degree;
                    cout << " ";
                }
            }
        }
        if (constant_term_sum >= 0) cout << "+";
        cout << constant_term_sum << " = 0" << endl;
    }

    // Generate initial seeds in the range [-20, 20] for both real and imaginary parts
    std::vector<std::complex<double>> initial_seeds = generateComplexInitialSeeds(num_vars, -20, 20);

    // Display the initial guesses
    cout << "\nInitial Guesses:\n";
    for (int i = 0; i < num_vars; ++i) 
    {
        char variable_name = 'a' + i;
        cout << "Variable " << variable_name << ": " << initial_seeds[i] << endl;
    }

    try {
        std::vector<std::complex<double>> solution = newtonsMethod(system, initial_seeds, 1e-2, 10000);
        cout << "\nSolution found:\n";
        for (int i = 0; i < num_vars; ++i) {
            cout << "Variable " << i + 1 << ": " << solution[i] << endl;
        }
    } catch (const std::exception& e) {
        cerr << e.what() << endl;
    }

    return 0;
}
