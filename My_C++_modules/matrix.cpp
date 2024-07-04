#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include <cmath>

using namespace std;

class Matrix {
private:
    vector<vector<double>> data;
    int rows, cols;

public:
    Matrix(int rows, int cols) : rows(rows), cols(cols) {
        data.resize(rows, vector<double>(cols, 0));
    }

    Matrix(const vector<vector<double>>& values) {
        rows = values.size();
        cols = values[0].size();
        data = values;
    }

    int getRows() const { return rows; }
    int getCols() const { return cols; }

    vector<double>& operator[](int i) { return data[i]; }
    const vector<double>& operator[](int i) const { return data[i]; }

    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result[j][i] = data[i][j];
        return result;
    }

    double determinant() const {
        if (rows != cols)
            throw invalid_argument("Matrix must be square to compute determinant.");

        Matrix temp(*this);
        double det = 1;
        for (int i = 0; i < rows; ++i) {
            int pivot = i;
            for (int j = i + 1; j < rows; ++j)
                if (abs(temp[j][i]) > abs(temp[pivot][i]))
                    pivot = j;

            if (temp[pivot][i] == 0)
                return 0;

            if (i != pivot) {
                swap(temp[i], temp[pivot]);
                det = -det;
            }

            det *= temp[i][i];
            for (int j = i + 1; j < rows; ++j)
                temp[i][j] /= temp[i][i];

            for (int j = 0; j < rows; ++j)
                if (j != i)
                    for (int k = i + 1; k < rows; ++k)
                        temp[j][k] -= temp[j][i] * temp[i][k];
        }
        return det;
    }

    Matrix inverse() const {
        if (rows != cols)
            throw invalid_argument("Matrix must be square to compute inverse.");

        int n = rows;
        Matrix result(n, n);
        Matrix temp(*this);

        for (int i = 0; i < n; ++i)
            result[i][i] = 1;

        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int j = i + 1; j < n; ++j)
                if (abs(temp[j][i]) > abs(temp[pivot][i]))
                    pivot = j;

            if (temp[pivot][i] == 0)
                throw runtime_error("Matrix is singular and cannot be inverted.");

            swap(temp[i], temp[pivot]);
            swap(result[i], result[pivot]);

            double pivotValue = temp[i][i];
            for (int j = 0; j < n; ++j) {
                temp[i][j] /= pivotValue;
                result[i][j] /= pivotValue;
            }

            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    double factor = temp[j][i];
                    for (int k = 0; k < n; ++k) {
                        temp[j][k] -= factor * temp[i][k];
                        result[j][k] -= factor * result[i][k];
                    }
                }
            }
        }
        return result;
    }

    static Matrix add(const Matrix& a, const Matrix& b) 
    {
        if (a.rows != b.rows || a.cols != b.cols)
            throw invalid_argument("Matrices must have the same dimensions to be added.");

        Matrix result(a.rows, a.cols);
        for (int i = 0; i < a.rows; ++i)
            for (int j = 0; j < a.cols; ++j)
                result[i][j] = a[i][j] + b[i][j];

        return result;
    }

    static Matrix subtract(const Matrix& A, const Matrix& B) 
    {
        if (A.rows != B.rows || A.cols != B.cols) {
            throw invalid_argument("Matrices must have the same dimensions for subtraction.");
        }

        Matrix result(A.rows, A.cols);

        for (int i = 0; i < A.rows; ++i) {
            for (int j = 0; j < A.cols; ++j) {
                result[i][j] = A[i][j] - B[i][j];
            }
        }

        return result;
    }

    static Matrix multiply(const Matrix& a, const Matrix& b) {
        if (a.cols != b.rows)
            throw invalid_argument("Number of columns of the first matrix must equal number of rows of the second matrix.");

        Matrix result(a.rows, b.cols);
        for (int i = 0; i < a.rows; ++i)
            for (int j = 0; j < b.cols; ++j)
                for (int k = 0; k < a.cols; ++k)
                    result[i][j] += a[i][k] * b[k][j];

        return result;
    }

    static Matrix identity(int size) 
    {
        Matrix id(size, size);
        for (int i = 0; i < size; ++i) {
            id[i][i] = 1.0;
        }
        return id;
    }

    static bool isDiagonal(const Matrix& mat) 
    {
        int size = mat.getRows();
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (i != j && mat[i][j] != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    static Matrix invertMatrix(const Matrix& A) 
    {
        int n = A.getRows();
        if (n != A.getCols()) {
            throw invalid_argument("Matrix must be square to invert.");
        }

        if (n == 1) {
            // Directly invert 1x1 matrix
            Matrix inv(1, 1);
            inv[0][0] = 1 / A[0][0];
            return inv;
        }

        // Determine block sizes
        int k = n - 1;
        Matrix A11(1, 1), A12(1, k), A21(k, 1), A22(k, k);

        // Extract blocks A11, A12, A21, A22
        A11[0][0] = A[0][0];
        for (int i = 0; i < k; ++i) {
            A12[0][i] = A[0][i + 1];
            A21[i][0] = A[i + 1][0];
            for (int j = 0; j < k; ++j) {
                A22[i][j] = A[i + 1][j + 1];
            }
        }

        // Invert A11
        Matrix A11_inv(1, 1);
        A11_inv[0][0] = 1 / A11[0][0];

        // Compute Schur complement S = A22 - A21 * A11_inv * A12
        Matrix temp = multiply(A21, multiply(A11_inv, A12));
        Matrix S = subtract(A22,temp);//   A22 - temp;

        // Invert S
        Matrix S_inv = invertMatrix(S);

        // Compute blocks of inverse matrix
        Matrix invA(n, n);
        // Upper-left block
        invA[0][0] = A11_inv[0][0];
        invA[0][0] = invA[0][0] + multiply(A11_inv,multiply(A12,multiply(S_inv,multiply(A21,A11_inv))))[0][0];

        // Upper-right block
        Matrix upper_right = multiply(A11_inv,multiply(A12,S_inv));
        for (int i = 0; i < k; ++i) 
        {
            invA[0][i + 1] = -upper_right[0][i];   //-A11_inv[0][0] * temp[0][i];
        }
        // Lower-left block
        Matrix lower_left = multiply(S_inv,multiply(A21,A11_inv));
        for (int i = 0; i < k; ++i) 
        {
            invA[i + 1][0] = -lower_left[i][0];   //-S_inv[i][0] * A11_inv[0][0];
        }
        // Lower-right block
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                invA[i + 1][j + 1] = S_inv[i][j];
            }
        }

        return invA;
    }


    void print() const {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                cout << setw(10) << data[i][j] << " ";
            }
            cout << endl;
        }
    }
};


int main() {
    vector<vector<double>> values;
    // = {
    //     {4, 7, 2},
    //     {3, 6, 1},
    //     {2, 5, 1}
    // };

    int m, n;

    cout << "Enter the matrix order m, n: ";
    cin >> m >> n;

    values.resize(m, vector<double>(n)); // Resize the values vector

    cout << "Enter the matrix elements row by row:" << endl;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> values[i][j]; // Use >> to read values
        }
    }

    Matrix mat(values);

    cout << "Original Matrix:" << endl;
    mat.print();

    cout << "Transpose:" << endl;
    Matrix transposed = mat.transpose();
    transposed.print();

    if (m == n) {
        cout << "Determinant: " << mat.determinant() << endl;

        try 
        {
            cout << "Inverse:" << endl;
            Matrix inversed =  Matrix::invertMatrix(mat);//  mat.inverse();
            inversed.print();

            Matrix inversed1 = mat.inverse();
            inversed1.print();
        } 
        catch (const runtime_error& e) {
            cout << e.what() << endl;
        }
    } 
    else 
    {
        cout << "Matrix is not square, so determinant and inverse cannot be computed." << endl;
    }

    return 0;
}