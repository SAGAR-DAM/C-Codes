#include <iostream>
#include <vector>
#include <complex>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include <utility>

using namespace std;

using Complex = complex<double>;

class Matrix {
private:
    vector<vector<Complex>> data;
    int rows, cols;

public:
    Matrix(int rows, int cols) : rows(rows), cols(cols) {
        data.resize(rows, vector<Complex>(cols, 0));
    }

    Matrix(const vector<vector<Complex>>& values) {
        rows = values.size();
        cols = values[0].size();
        data = values;
    }

    // New constructor to read matrix from cin
    Matrix() {
        cout << "Enter number of rows: ";
        cin >> rows;
        cout << "Enter number of columns: ";
        cin >> cols;
        data.resize(rows, vector<Complex>(cols, 0));
        cout << "Enter matrix elements row by row (real and imaginary parts):" << endl;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double real, imag;
                cin >> real >> imag;
                data[i][j] = Complex(real, imag);
            }
        }
    }

    int getRows() const { return rows; }
    int getCols() const { return cols; }

    vector<Complex>& operator[](int i) { return data[i]; }
    const vector<Complex>& operator[](int i) const { return data[i]; }

    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result[j][i] = data[i][j];
        return result;
    }

    Complex determinant() const {
        if (rows != cols)
            throw invalid_argument("Matrix must be square to compute determinant.");

        Matrix temp(*this);
        Complex det = 1;
        for (int i = 0; i < rows; ++i) {
            int pivot = i;
            for (int j = i + 1; j < rows; ++j)
                if (abs(temp[j][i]) > abs(temp[pivot][i]))
                    pivot = j;

            if (temp[pivot][i] == Complex(0, 0))
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

            if (temp[pivot][i] == Complex(0, 0))
                throw runtime_error("Matrix is singular and cannot be inverted.");

            swap(temp[i], temp[pivot]);
            swap(result[i], result[pivot]);

            Complex pivotValue = temp[i][i];
            for (int j = 0; j < n; ++j) {
                temp[i][j] /= pivotValue;
                result[i][j] /= pivotValue;
            }

            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    Complex factor = temp[j][i];
                    for (int k = 0; k < n; ++k) {
                        temp[j][k] -= factor * temp[i][k];
                        result[j][k] -= factor * result[i][k];
                    }
                }
            }
        }
        return result;
    }

    static Matrix add(const Matrix& a, const Matrix& b) {
        if (a.rows != b.rows || a.cols != b.cols)
            throw invalid_argument("Matrices must have the same dimensions to be added.");

        Matrix result(a.rows, a.cols);
        for (int i = 0; i < a.rows; ++i)
            for (int j = 0; j < a.cols; ++j)
                result[i][j] = a[i][j] + b[i][j];

        return result;
    }

    static Matrix subtract(const Matrix& A, const Matrix& B) {
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

    static Matrix identity(int size) {
        Matrix id(size, size);
        for (int i = 0; i < size; ++i) {
            id[i][i] = 1.0;
        }
        return id;
    }

    static bool isDiagonal(const Matrix& mat) {
        int size = mat.getRows();
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (i != j && mat[i][j] != Complex(0, 0)) {
                    return false;
                }
            }
        }
        return true;
    }

    static Matrix invertMatrix(Matrix A) {
        int n = A.getRows();
        if (n != A.getCols()) {
            throw invalid_argument("Matrix must be square to invert.");
        }

        if (n == 1) {
            // Directly invert 1x1 matrix
            Matrix inv(1, 1);
            inv[0][0] = Complex(1,0) / A[0][0];
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

        Matrix A11_inv = A11.inverse();
        Matrix S = subtract(A22, multiply(A21, multiply(A11_inv, A12)));

        Matrix S_inv = S.inverse();

        Matrix B11 = add(A11_inv, multiply(multiply(A11_inv, A12), multiply(S_inv, A21)));
        Matrix B12 = multiply(multiply(A11_inv, A12), S_inv);
        Matrix B21 = multiply(multiply(S_inv, A21), A11_inv);
        Matrix B22 = S_inv;

        Matrix result(n, n);
        result[0][0] = B11[0][0];
        for (int i = 0; i < k; ++i) {
            result[0][i + 1] = B12[0][i];
            result[i + 1][0] = B21[i][0];
            for (int j = 0; j < k; ++j) {
                result[i + 1][j + 1] = Complex(-1,0)*B22[i][j];
            }
        }

        return result;
    }

    static pair<Matrix, Matrix> qrDecomposition(Matrix A) 
    {
        int n = A.getRows();
        Matrix Q = Matrix::identity(n); // Initialize Q to identity matrix
        Matrix R = A;

        for (int k = 0; k < n - 1; ++k) {
            Matrix H = Matrix::identity(n);
            Complex norm_x = 0;
            for (int i = k; i < n; ++i) {
                norm_x += norm(R[i][k]);
            }
            norm_x = sqrt(norm_x);

            Complex rkk = R[k][k];
            Complex s = (real(rkk) >= 0) ? 1.0 : -1.0;
            Complex u1 = rkk + s * norm_x;
            Complex w = sqrt(2.0 * norm_x * (norm_x + abs(rkk)));

            for (int i = k; i < n; ++i) {
                H[i][k] = R[i][k] / w;
            }
            H[k][k] = 1;

            R = multiply(H, R);
            Q = multiply(Q, H.transpose());
        }

        return make_pair(Q, R);
    }

static vector<Complex> eigenvalues(Matrix A, int maxIterations = 1000, double tolerance = 1e-10) 
    {
        int n = A.getRows();
        
        if (n != A.getCols()) {
            throw invalid_argument("Matrix must be square.");
        }

        Matrix Ak = A;
        Matrix Q = Matrix::identity(n), R = Matrix::identity(n);

        for (int i = 0; i < maxIterations; ++i) {
            tie(Q, R) = qrDecomposition(Ak);
            Ak = multiply(R, Q);

            // Check convergence
            bool converged = true;
            for (int j = 0; j < n - 1; ++j) {
                if (abs(Ak[j + 1][j]) > tolerance) {
                    converged = false;
                    break;
                }
            }
            if (converged) break;
        }

        vector<Complex> eigenvalues(n);
        for (int i = 0; i < n; ++i) {
            eigenvalues[i] = Ak[i][i];
        }

        return eigenvalues;
    }

};

int main() {
    vector<vector<Complex>> values 
    ={
        {Complex(1,0), Complex(-1,0)},
        {Complex(1,0), Complex(1,0)}
    };

    // = {
    //     {Complex(4, 1), Complex(2, 7), Complex(0, 0)},
    //     {Complex(4, 9), Complex(6, 4), Complex(0, 0)},
    //     {Complex(0, 0), Complex(0, 0), Complex(1, 7)}
    // };


    // Uncomment below lines for user input
    // int m, n;
    // cout << "Enter the matrix order m, n: ";
    // cin >> m >> n;
    // values.resize(m, vector<Complex>(n)); // Resize the values vector
    // cout << "Enter the matrix elements row by row (real and imaginary parts):" << endl;
    // for (int i = 0; i < m; ++i) {
    //     for (int j = 0; j < n; ++j) {
    //         double real, imag;
    //         cin >> real >> imag;
    //         values[i][j] = Complex(real, imag);
    //     }
    // }

    Matrix mat(values);

    cout << "Matrix elements:" << endl;
    for (int i = 0; i < mat.getRows(); ++i) {
        for (int j = 0; j < mat.getCols(); ++j) {
            cout << "(" << mat[i][j].real() << ", " << mat[i][j].imag() << ") ";
        }
        cout << endl;
    }

    // try {
    //     Complex det = mat.determinant();
    //     cout << "Determinant: " << det << endl;
    //     Matrix inv = mat.inverse();
    //     cout << "Inverse matrix elements:" << endl;
    //     for (int i = 0; i < inv.getRows(); ++i) {
    //         for (int j = 0; j < inv.getCols(); ++j) {
    //             cout << "(" << inv[i][j].real() << ", " << inv[i][j].imag() << ") ";
    //         }
    //         cout << endl;
    //     }
    // } catch (const exception& e) {
    //     cerr << e.what() << endl;
    // }


    vector<Complex> eig = Matrix::eigenvalues(mat);
    for(int i=0; i<eig.size(); i++)
    {
        cout<<eig[i]<<", ";
    }

    return 0;
}
