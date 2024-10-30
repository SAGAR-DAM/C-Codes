#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <chrono>

using namespace std;

// Function to generate a random integer matrix (2D vector)
vector<vector<double>> generateRandomDoubleMatrix(int rows, int cols, double min_value, double max_value) {
    unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
    mt19937 gen(seed); // Mersenne Twister engine
    uniform_real_distribution<> distrib(min_value, max_value);

    vector<vector<double>> random_matrix(rows, vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            random_matrix[i][j] = distrib(gen); // Generate a random value for each element
        }
    }
    return random_matrix;
}


// Function to generate a random integer matrix (2D vector) with double type entries
vector<vector<double>> generateRandomIntMatrix(int rows, int cols, int min_value, int max_value) {
    unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
    mt19937 gen(seed); // Mersenne Twister engine
    uniform_int_distribution<int> distrib(min_value, max_value);

    vector<vector<double>> random_matrix(rows, vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            random_matrix[i][j] = static_cast<double>(distrib(gen)); // Generate a random integer and cast to double
        }
    }
    return random_matrix;
}


// Function to perform Gaussian elimination to solve Ax = B
vector<double> gaussianElimination(vector<vector<double>>& A, vector<double>& B) {
    int n = A.size();
    
    // Forward elimination
    for (int i = 0; i < n; i++) {
        // Partial pivoting
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        // Swap rows in A and B
        swap(A[i], A[maxRow]);
        swap(B[i], B[maxRow]);

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; k++) {
            double c = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= c * A[i][j];
            }
            B[k] -= c * B[i];
        }
    }

    // Back substitution
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = B[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}

// Function to print the system of equations
void printSystem(vector<vector<double>>& A, vector<double>& B) {
    int n = A.size();
    cout << "System of equations:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (A[i][j] >= 0 && j > 0) cout << "+";
            cout << A[i][j] << "x" << (j + 1) << " ";
        }
        cout << "= " << B[i] << endl;
    }
    cout << endl;
}

int main() {
    int num_variables;
    cout << "Enter the number of variables: ";
    cin >> num_variables;

    // Generate random matrix A and vector B
    vector<vector<double>> A = generateRandomIntMatrix(num_variables, num_variables, -10, 10);
    vector<double> B(num_variables);
    for (int i = 0; i < num_variables; i++) {
        B[i] = generateRandomIntMatrix(1, 1, -10 * num_variables, 10 * num_variables)[0][0];
    }

    cout<<"\033[31m";
    cout<<"================================================================="<<endl;
    printSystem(A, B);
    cout<<"================================================================="<<endl;

    vector<double> solution = gaussianElimination(A, B);

    cout << "\033[34m" << "Solution:\n------------\n"<<"\033[32m";
    for (int i = 0; i < num_variables; i++) {
        cout << "x" << (i + 1) << " = " << solution[i] << endl;
    }
    cout <<"================================================================="<<endl;
    cout << "\033[0m" << "Press any key to exit" << endl;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cin.get();
    return 0;
}
