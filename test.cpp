#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

// Function to perform the Gram-Schmidt process
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>> gramSchmidt(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    std::vector<std::vector<double>> Q(n, std::vector<double>(n, 0));
    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0));

    for (int j = 0; j < n; ++j) {
        std::vector<double> v = A[j];
        for (int i = 0; i < j; ++i) {
            R[i][j] = 0;
            for (int k = 0; k < n; ++k) {
                R[i][j] += Q[k][i] * A[k][j];
            }
            for (int k = 0; k < n; ++k) {
                v[k] -= R[i][j] * Q[k][i];
            }
        }
        R[j][j] = 0;
        for (int k = 0; k < n; ++k) {
            R[j][j] += v[k] * v[k];
        }
        R[j][j] = std::sqrt(R[j][j]);
        for (int k = 0; k < n; ++k) {
            Q[k][j] = v[k] / R[j][j];
        }
    }

    return std::make_tuple(Q, R);
}

// Function to multiply two matrices
std::vector<std::vector<double>> matMul(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    int n = A.size();
    std::vector<std::vector<double>> C(n, std::vector<double>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// QR Algorithm to find all eigenvalues
std::vector<double> qrAlgorithm(std::vector<std::vector<double>> A, int maxIterations = 1000, double tolerance = 1e-10) {
    int n = A.size();
    for (int iter = 0; iter < maxIterations; ++iter) {
        auto [Q, R] = gramSchmidt(A);
        A = matMul(R, Q);
    }

    std::vector<double> eigenvalues(n);
    for (int i = 0; i < n; ++i) {
        eigenvalues[i] = A[i][i];
    }

    return eigenvalues;
}

int main() {
    // Example usage
    std::vector<std::vector<double>> matrix = {
        {4, 1, 2},
        {1, 3, 0},
        {2, 0, 1}
    };

    std::vector<double> eigenvalues = qrAlgorithm(matrix);

    std::cout << "Eigenvalues: ";
    for (double eigenvalue : eigenvalues) {
        std::cout << eigenvalue << " ";
    }
    std::cout << std::endl;

    return 0;
}
