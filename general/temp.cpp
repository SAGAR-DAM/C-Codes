#include<iostream>
#include<string>
#include<limits>
#include<vector>
#include<algorithm>
#include "D:\C++ codes\My_C++_modules\np.hpp"
#include <stdexcept> // For std::invalid_argument

using namespace std;

class Matrix {
private:
    std::vector<std::vector<double>> data;
    size_t rows;
    size_t cols;

public:
    // Constructors
    Matrix(size_t rows, size_t cols, double initial = 0.0)
        : rows(rows), cols(cols), data(std::vector<std::vector<double>>(rows, std::vector<double>(cols, initial))) {}

    // Accessors
    size_t numRows() const { return rows; }
    size_t numCols() const { return cols; }

    // Element access
    double& operator()(size_t row, size_t col) {
        if (row >= rows || col >= cols) {
            throw std::out_of_range("Matrix indices out of range");
        }
        return data[row][col];
    }

    double operator()(size_t row, size_t col) const {
        if (row >= rows || col >= cols) {
            throw std::out_of_range("Matrix indices out of range");
        }
        return data[row][col];
    }

    // Matrix addition
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for addition");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] + other(i, j);
            }
        }
        return result;
    }

    // Matrix subtraction
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("Matrix dimensions must match for subtraction");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] - other(i, j);
            }
        }
        return result;
    }

    // Scalar multiplication
    Matrix operator*(double scalar) const {
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] * scalar;
            }
        }
        return result;
    }

    // Scalar division
    Matrix operator/(double scalar) const {
        if (scalar == 0.0) {
            throw std::invalid_argument("Division by zero");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] / scalar;
            }
        }
        return result;
    }

    // Matrix multiplication
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Matrix dimensions are incompatible for multiplication");
        }
        Matrix result(rows, other.cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.cols; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < cols; ++k) {
                    sum += data[i][k] * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    // Matrix exponentiation (only for square matrices and positive integer exponents)
    Matrix operator^(int exp) const {
        if (rows != cols) {
            throw std::invalid_argument("Matrix must be square for exponentiation");
        }
        if (exp < 0) {
            throw std::invalid_argument("Exponent must be non-negative");
        }
        Matrix result(rows, cols);
        Matrix base = *this;
        result.identity(); // Initialize result as identity matrix
        while (exp > 0) {
            if (exp % 2 == 1) {
                result = result * base;
            }
            base = base * base;
            exp /= 2;
        }
        return result;
    }

    // Helper function to initialize as identity matrix (for exponentiation)
    void identity() {
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                if (i == j) {
                    data[i][j] = 1.0;
                } else {
                    data[i][j] = 0.0;
                }
            }
        }
    }

    // Output operator
    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
        for (size_t i = 0; i < matrix.rows; ++i) {
            for (size_t j = 0; j < matrix.cols; ++j) {
                os << matrix(i, j) << " ";
            }
            os << "\n";
        }
        return os;
    }
};

int main() {
    // Example usage
    Matrix A(2, 2, 1.0);
    Matrix B(2, 2, 2.0);

    std::cout << "Matrix A:\n" << A << std::endl;
    std::cout << "Matrix B:\n" << B << std::endl;

    Matrix C = A + B;
    Matrix D = A - B;
    Matrix E = A * 2.0;
    Matrix F = A / 2.0;
    Matrix G = A * B;

    std::cout << "Matrix C (A + B):\n" << C << std::endl;
    std::cout << "Matrix D (A - B):\n" << D << std::endl;
    std::cout << "Matrix E (A * 2.0):\n" << E << std::endl;
    std::cout << "Matrix F (A / 2.0):\n" << F << std::endl;
    std::cout << "Matrix G (A * B):\n" << G << std::endl;

    // Example of matrix exponentiation (A^2)
    Matrix H(2, 2);
    H(0, 0) = 1.0;
    H(0, 1) = 2.0;
    H(1, 0) = 3.0;
    H(1, 1) = 4.0;
    Matrix I = H ^ 2;
    std::cout << "Matrix I (H^2):\n" << I << std::endl;

    return 0;
}
