#include<iostream>
#include<string>
#include<limits>
#include"D:\C++ codes\My_C++_modules\np.hpp"
#include<cmath>
#include<algorithm>
#include<complex>

using namespace std;
using namespace std;

int main() {
    // Read matrix dimensions from input (assuming it's passed by Python)
    int rows, cols;
    cin >> rows >> cols;

    // Read matrix elements
    vector<vector<double>> values(rows, vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cin >> values[i][j];
        }
    }

    // Create a Matrix object
    Matrix mat(values);

    // Compute coefficients of characteristic polynomial
    vector<double> coefficients = Matrix::coeff_charac(mat);

    // Print coefficients to standard output
    for (int i = 0; i < coefficients.size(); ++i) {
        cout << coefficients[i] << " ";
    }
    cout << endl;

    return 0;
}
