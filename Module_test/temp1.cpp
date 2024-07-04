#include <iostream>
#include <vector>
#include <iomanip> // For std::setw
#include <cmath>

using namespace std;

// Function to perform the Gaussian Elimination
bool gaussianElimination(vector<vector<double>>& matrix, vector<vector<double>>& inverse) {
    int n = matrix.size();

    // Augment the matrix with the identity matrix
    for (int i = 0; i < n; ++i) {
        inverse[i][i] = 1;
    }

    for (int i = 0; i < n; ++i) {
        // Search for maximum in this column
        double maxEl = abs(matrix[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(matrix[k][i]) > maxEl) {
                maxEl = abs(matrix[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k = 0; k < n; ++k) {
            swap(matrix[maxRow][k], matrix[i][k]);
            swap(inverse[maxRow][k], inverse[i][k]);
        }

        // Make sure the matrix is not singular
        if (matrix[i][i] == 0) {
            return false; // Singular matrix
        }

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; ++k) {
            double c = -matrix[k][i] / matrix[i][i];
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    matrix[k][j] = 0;
                } else {
                    matrix[k][j] += c * matrix[i][j];
                }
                inverse[k][j] += c * inverse[i][j];
            }
        }
    }

    // Solve equation Ax = b for an upper triangular matrix A
    for (int i = n - 1; i >= 0; --i) {
        for (int j = 0; j < n; ++j) {
            inverse[i][j] /= matrix[i][i];
        }
        for (int k = i - 1; k >= 0; --k) {
            double c = -matrix[k][i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    matrix[k][j] += c * matrix[i][j];
                }
                inverse[k][j] += c * inverse[i][j];
            }
        }
    }
    return true;
}

// Function to print a matrix
void printMatrix(const vector<vector<double>>& matrix) {
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << setw(10) << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

int main() {
    int n;
    cout << "Enter the size of the matrix (n x n): ";
    cin >> n;

    vector<vector<double>> matrix(n, vector<double>(n));
    vector<vector<double>> inverse(n, vector<double>(n, 0));

    cout << "Enter the elements of the matrix row by row:" << endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> matrix[i][j];
        }
    }

    cout << "Original Matrix:" << endl;
    printMatrix(matrix);

    if (gaussianElimination(matrix, inverse)) {
        cout << "Inverse Matrix:" << endl;
        printMatrix(inverse);
    } else {
        cout << "Matrix is singular and cannot be inverted." << endl;
    }

    return 0;
}
