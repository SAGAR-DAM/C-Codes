#include<iostream>
#include<string>
#include<limits>
#include<vector>
#include<algorithm>
#include<iomanip>
#include <stdexcept>
#include "D:\C++ codes\My_C++_modules\np.hpp"
#include<cmath>

using namespace std;

/*-------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------- */
//                          FUNCTION DECLARATION:
/*-------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------- */



// ###########################################################################
//          Function to find the index of a value in a vector
// ###########################################################################
int find_index(const std::vector<double>& array, double value) {
    // Iterate through the vector to find the value
    for (int i = 0; i < array.size(); ++i) {
        if (array[i] == value) {
            return i; // Return the index if found
        }
    }
    // Return -1 if the value is not found in the vector
    cout<<"Item not found...\n";
    return -1;
}


// ###########################################################################
//  Function to find the index of the closest value in a vector of doubles
// ###########################################################################
int find_closest_index(const std::vector<double>& array, double value) {
    // Check if the vector is empty
    if (array.empty()) {
        return -1; // Return -1 if the vector is empty
    }

    // Initialize variables to track the closest index and the minimum difference
    int closestIndex = 0; // Assume the closest index is the first element initially
    double minDifference = std::abs(array[0] - value); // Initialize with difference to first element

    // Iterate through the vector to find the index of the closest value
    for (int i = 1; i < array.size(); ++i) {
        double difference = std::abs(array[i] - value);
        if (difference < minDifference) {
            minDifference = difference;
            closestIndex = i;
        }
    }
    //cout<<array[closestIndex]<<endl;
    return closestIndex;
}


// ###########################################################################
//                      cut array for given range
// ###########################################################################
std::vector<double> cutarr(std::vector<double> array, int minw, int maxw)
{
    vector<double> newarray(array.begin() + minw, array.begin()+maxw);
    return newarray;
}

// ###########################################################################
//           Function to calculate the mean of a vector of doubles
// ###########################################################################
double meanarr(const std::vector<double>& vec) {
    // Check if the vector is empty to avoid division by zero
    if (vec.empty()) {
        return 0.0; // Return 0 if the vector is empty
    }

    // Calculate the sum of all elements in the vector
    double sum = 0.0;
    for (double num : vec) {
        sum += num;
    }

    // Calculate the mean (average)
    double mean = sum / vec.size();
    return mean;
}

// ###########################################################################
//      Function to calculate the standard deviation of a vector of doubles
// ###########################################################################
double stdarr(const std::vector<double>& vec) {
    if (vec.empty()) {
        return 0.0; // Return 0 if the vector is empty
    }

    double mean = meanarr(vec); // Calculate mean using the previously defined function

    // Calculate sum of squared differences from the mean
    double sumSquaredDiff = 0.0;
    for (double num : vec) {
        sumSquaredDiff += (num - mean) * (num - mean);
    }

    // Calculate variance
    double variance = sumSquaredDiff / vec.size();

    // Calculate standard deviation (square root of variance)
    double stdDev = std::sqrt(variance);
    return stdDev;
}




/*-------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------- */
//                             CLASS DECLARATION:
/*-------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------- */
/*-------------------------------------------------------------------------- */


// ###########################################################################
//                              COMPLEX CLASS
// ###########################################################################

// costructor of the class: 
//----------------------------
complex::complex(double r, double i) : real(r), imag(i) {}


// Functions over the class:
//-----------------------------
double complex::getreal()
{
    return real;
}

double complex::getimag()
{
    return imag;
}

double complex::get_r()
{
    return r;
}

float complex::get_theta()
{
    return theta;
}

void complex::getcomplex()
{
    if (imag >= 0)
    {
        cout << "The complex number is: " << real << " + " << imag << "j" << endl;
    }
    else
    {
        cout << "The complex number is: " << real << " " << imag << "j" << endl;
    }
}

void complex::getcomplex_polar()
{
    //rect_to_polar();
    cout << "The complex number in polar form: " << r << " exp(" << theta << "j)" << endl;
}

double complex::abscomplex(complex c)
{
    return pow((c.getreal()*c.getreal()+c.getimag()*c.getimag()),0.5);
}

float complex::argcomplex(complex c)
{   
    //c.rect_to_polar();
    return c.theta; 
}

complex complex::addcomplex(complex c1, complex c2)
{
    return complex(c1.real + c2.real, c1.imag + c2.imag);
}

complex complex::subtractcomplex(complex c1, complex c2)
{
    return complex(c1.real - c2.real, c1.imag - c2.imag);
}

complex complex::multcomplex(complex c1, complex c2)
{
    //c1.rect_to_polar();
    //c2.rect_to_polar();
    double r_mult = c1.r * c2.r;
    float theta_mult = (c1.theta + c2.theta);

    if (theta_mult > M_PI)
    {
        theta_mult = -(2 * M_PI - theta_mult);
    }
    if (theta_mult < M_PI)
    {
        theta_mult = 2 * M_PI + theta_mult;
    }

    double real_mult = r_mult * cos(theta_mult);
    double imag_mult = r_mult * sin(theta_mult);

    return complex(real_mult, imag_mult);
}

complex complex::divcomplex(complex c1, complex c2)
{
    //c1.rect_to_polar();
    //c2.rect_to_polar();
    double r_div = c1.r / c2.r;
    float theta_div = (c1.theta - c2.theta);

    if (theta_div > M_PI)
    {
        theta_div = -(2 * M_PI - theta_div);
    }
    if (theta_div < M_PI)
    {
        theta_div = 2 * M_PI + theta_div;
    }

    double real_div = r_div * cos(theta_div);
    double imag_div = r_div * sin(theta_div);

    return complex(real_div, imag_div);
}

complex complex::complex_exponenciation(complex c, float n)
{
    //c.rect_to_polar();
    double r_exp = pow(c.r, n);
    float theta_exp = n * c.theta;

    double real_exp = r_exp * cos(theta_exp);
    double imag_exp = r_exp * sin(theta_exp);
    return complex(real_exp, imag_exp);
}

complex complex::complexpow(complex c1, complex c2)
{
    //c1.rect_to_polar();
    //c2.rect_to_polar();
    double r1 = c1.r, r2 = c2.r;
    float theta1 = c1.theta, theta2 = c2.theta;

    double r_pow = (pow(r1, r2 * cos(theta2))) * exp(-r2 * theta1 * sin(theta2));
    float theta_pow = r2 * sin(theta2) * log(r1) + r2 * cos(theta2) * theta1;

    double real_pow = r_pow * cos(theta_pow);
    double imag_pow = r_pow * sin(theta_pow);
    return complex(real_pow, imag_pow);
}

complex complex::real_pow_complex(complex c, float n)
{
    complex c1 = complexpow(complex(n,0),c);
    return c1;
}

complex complex::sincomplex(complex c)
{
    double real_sin = sin(c.real)*cosh(c.imag), imag_sin = cos(c.real)*sinh(c.imag);
    complex c_sin = complex(real_sin, imag_sin);
    //c_sin.rect_to_polar();
    return c_sin;
}

complex complex::coscomplex(complex c)
{
    double real_cos = cos(c.real)*cosh(c.imag), imag_cos = -sin(c.real)*sinh(c.imag);
    complex c_cos = complex(real_cos, imag_cos);
    //c_cos.rect_to_polar();
    return c_cos;
}

complex complex::tancomplex(complex c)
{
    complex c_tan = divcomplex(sincomplex(c),coscomplex(c));
    return c_tan;
}

complex complex::csccomplex(complex c)
{
    complex c_csc = divcomplex(complex(1,0),sincomplex(c));
    return c_csc;
}

complex complex::seccomplex(complex c)
{
    complex c_sec = divcomplex(complex(1,0),coscomplex(c));
    return c_sec;
}

complex complex::cotcomplex(complex c)
{
    complex c_cot = divcomplex(coscomplex(c),sincomplex(c));
    return c_cot;
}

complex complex::logcomplex(complex c)
{
    //c.rect_to_polar();
    double real_log = log(c.r), imag_log = c.theta;
    complex c_log = complex(real_log,imag_log);
    return c_log;
}

complex complex::asincomplex(complex c)
{
    complex iota = complex(0,1);
    complex z = addcomplex(multcomplex(iota,c),complex_exponenciation(subtractcomplex(complex(1,0),multcomplex(c,c)),0.5));
    z = logcomplex(z);
    z = divcomplex(z,iota);
    return z;
}

complex complex::acoscomplex(complex c)
{
    complex iota = complex(0,1);
    complex z = addcomplex(c,complex_exponenciation(subtractcomplex(multcomplex(c,c),complex(1,0)),0.5));
    z = logcomplex(z);
    z = divcomplex(z,iota);
    return z;
}

complex complex::atancomplex(complex c)
{
    complex iota = complex(0,1);
    complex z = divcomplex(subtractcomplex(iota,c),addcomplex(iota,c));
    z = logcomplex(z);
    z = divcomplex(z,complex(0,2));
    return z;
}

complex complex::acsccomplex(complex c)
{
    complex z = asincomplex(divcomplex(complex(1,0),c));
    return z;
}

complex complex::aseccomplex(complex c)
{
    complex z = acoscomplex(divcomplex(complex(1,0),c));
    return z;
}

complex complex::acotcomplex(complex c)
{
    complex z = atancomplex(divcomplex(complex(1,0),c));
    return z;
}


// Addition operator definition outside the class
complex operator+(complex& c1, complex& c2) 
{
    return complex::addcomplex(c1,c2);
}


// Operators of the complex class:
//----------------------------------
complex operator-(complex& c1, complex& c2) 
{
    return complex::subtractcomplex(c1,c2);
}

// Multiplication operator definition outside the class
complex operator*(complex& c1, complex& c2) 
{
    return complex::multcomplex(c1,c2);
}

// Division operator definition outside the class
complex operator/(complex& c1, complex& c2) 
{
    return complex::divcomplex(c1,c2);
}

// Exponentiation operator definition outside the class
complex operator^(complex& c1, complex& c2) 
{
    return complex::complexpow(c1,c2);
}

// Output operator
std::ostream& operator<<(std::ostream& os, complex& c) {
    os << c.getreal();
    if (c.getimag() >= 0) {
        os << " + " << c.getimag() << "j";
    } else {
        os << " - " << -c.getimag() << "j";
    }
    return os;
}




// ###########################################################################
//                              COMPLEX CLASS
// ###########################################################################

Matrix::Matrix(int rows, int cols) : rows(rows), cols(cols) 
{
    data.resize(rows, vector<double>(cols, 0));
}

Matrix::Matrix(vector<vector<double>>& values) 
{
    rows = values.size();
    cols = values[0].size();
    data = values;
}

int Matrix::getRows() { return rows; }
int Matrix::getCols() { return cols; }

vector<double>& Matrix::operator[](int i) { return data[i]; }
//const vector<double>& Matrix::operator[](int i) { return data[i]; }

Matrix Matrix::transpose() 
{
    Matrix result(cols, rows);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result[j][i] = data[i][j];
    return result;
}

double Matrix::determinant() 
{
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

Matrix Matrix::inverse()
{
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

Matrix Matrix::add(Matrix a, Matrix b) 
{
    if (a.rows != b.rows || a.cols != b.cols)
    {
        throw invalid_argument("Matrices must have the same dimensions to be added.");
    }

    Matrix result(a.rows, a.cols);
    for (int i = 0; i < a.rows; ++i)
        for (int j = 0; j < a.cols; ++j)
            result[i][j] = a[i][j] + b[i][j];

    return result;
}

Matrix Matrix::subtract(Matrix A, Matrix B) 
{
    if (A.rows != B.rows || A.cols != B.cols) 
    {
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

Matrix Matrix::multiply(Matrix a, Matrix b) 
{
    if (a.cols != b.rows)
    {
        throw invalid_argument("Number of columns of the first matrix must equal number of rows of the second matrix.");
    }
    Matrix result(a.rows, b.cols);
    for (int i = 0; i < a.rows; ++i)
        for (int j = 0; j < b.cols; ++j)
            for (int k = 0; k < a.cols; ++k)
                result[i][j] += a[i][k] * b[k][j];

    return result;
}

Matrix Matrix::identity(int size) 
{
    Matrix id(size, size);
    for (int i = 0; i < size; ++i) 
    {
        id[i][i] = 1.0;
    }
    return id;
}

bool Matrix::isDiagonal(Matrix mat) 
{
    int size = mat.getRows();
    for (int i = 0; i < size; ++i) 
    {
        for (int j = 0; j < size; ++j) 
        {
            if (i != j && mat[i][j] != 0) 
            {
                return false;
            }
        }
    }
    return true;
}


Matrix Matrix::invertMatrix(Matrix A) 
{
    int n = A.getRows();
    if (n != A.getCols()) 
    {
        throw invalid_argument("Matrix must be square to invert.");
    }

    if (n == 1) 
    {
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
    for (int i = 0; i < k; ++i) 
    {
        A12[0][i] = A[0][i + 1];
        A21[i][0] = A[i + 1][0];
        for (int j = 0; j < k; ++j) 
        {
            A22[i][j] = A[i + 1][j + 1];
        }
    }

    // Invert A11
    Matrix A11_inv(1, 1);
    A11_inv[0][0] = 1 / A11[0][0];

    // Compute Schur complement S = A22 - A21 * A11_inv * A12
    Matrix temp = Matrix::multiply(A21, Matrix::multiply(A11_inv, A12));
    Matrix S = Matrix::subtract(A22,temp);//   A22 - temp;

    // Invert S
    Matrix S_inv = invertMatrix(S);

    // Compute blocks of inverse matrix
    Matrix invA(n, n);
    // Upper-left block
    invA[0][0] = A11_inv[0][0];
    invA[0][0] = invA[0][0] + Matrix::multiply(A11_inv, Matrix::multiply(A12,Matrix::multiply(S_inv,Matrix::multiply(A21,A11_inv))))[0][0];

    // Upper-right block
    Matrix upper_right = Matrix::multiply(A11_inv,Matrix::multiply(A12,S_inv));
    for (int i = 0; i < k; ++i) 
    {
        invA[0][i + 1] = -upper_right[0][i];   //-A11_inv[0][0] * temp[0][i];
    }
    // Lower-left block
    Matrix lower_left = Matrix::multiply(S_inv,Matrix::multiply(A21,A11_inv));
    for (int i = 0; i < k; ++i) 
    {
        invA[i + 1][0] = -lower_left[i][0];   //-S_inv[i][0] * A11_inv[0][0];
    }
    // Lower-right block
    for (int i = 0; i < k; ++i) 
    {
        for (int j = 0; j < k; ++j) 
        {
            invA[i + 1][j + 1] = S_inv[i][j];
        }
    }

    return invA;
}


pair<Matrix, Matrix> Matrix::qrDecomposition(Matrix A) 
{
    int n = A.getRows();
    Matrix Q = Matrix::identity(n); // Initialize Q to identity matrix
    Matrix R = A;

    for (int k = 0; k < n - 1; ++k) 
    {
        Matrix H = Matrix::identity(n);
        double norm_x = 0;
        for (int i = k; i < n; ++i) 
        {
            norm_x += R[i][k] * R[i][k];
        }
        norm_x = sqrt(norm_x);

        double rkk = R[k][k];
        double s = (rkk >= 0) ? 1 : -1;
        double u1 = rkk + s * norm_x;
        double w = sqrt(2 * norm_x * (norm_x + abs(rkk)));
        
        for (int i = k; i < n; ++i) 
        {
            H[i][k] = R[i][k] / w;
        }
        H[k][k] = 1;

        R = Matrix::multiply(H , R);
        Q = Matrix::multiply(Q , H.transpose());
    }

    return make_pair(Q, R);
}

vector<double> Matrix::eigenvalues(Matrix A, int maxIterations = 1000, double tolerance = 1e-10) 
{
    int n = A.getRows();
    
    if (n != A.getCols()) 
    {
        throw invalid_argument("Matrix must be square.");
    }

    Matrix Ak = A;
    Matrix Q = Matrix::identity(n), R = Matrix::identity(n);

    for (int i = 0; i < maxIterations; ++i) 
    {
        tie(Q, R) = Matrix::qrDecomposition(Ak);
        Ak = Matrix::multiply(R, Q);

        // Check convergence
        bool converged = true;
        for (int j = 0; j < n - 1; ++j) 
        {
            if (abs(Ak[j + 1][j]) > tolerance) 
            {
                converged = false;
                break;
            }
        }
        if (converged) break;
    }

    vector<double> eigenvalues(n);
    for (int i = 0; i < n; ++i) 
    {
        eigenvalues[i] = Ak[i][i];
    }

    return eigenvalues;
}

Matrix Matrix::qrAlgorithm(Matrix A, int maxIterations = 1000, double tolerance = 1e-10) 
{
    int n = A.getRows();

    Matrix Q_total = Matrix::identity(n);
    Matrix Ak = A;
    Matrix Q = Matrix::identity(n), R = Matrix::identity(n);

    for (int i = 0; i < maxIterations; ++i) {
        tie(Q, R) = Matrix::qrDecomposition(Ak);
        Ak = Matrix::multiply(R, Q);
        Q_total = Matrix::multiply(Q_total, Q);

        // Check convergence
        bool converged = true;
        for (int j = 0; j < n; ++j) {
            if (abs(Ak[j][j] - Q_total[j][j]) > tolerance) {
                converged = false;
                break;
            }
        }
        if (converged) break;
    }

    return Q_total;
}

Matrix Matrix::similarity_transform(Matrix A) 
{
    if (A.getRows() != A.getCols()) {
        throw invalid_argument("Matrix must be square.");
    }

    // Perform QR algorithm to get the matrix of eigenvectors
    Matrix eigenvectorMatrix = qrAlgorithm(A);

    return eigenvectorMatrix;
}

void Matrix::print()
{
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << setw(10) << data[i][j] << " ";
        }
        std::cout << endl;
    }
}


// Operators of the complex class:
//----------------------------------
Matrix operator+(Matrix& M1, Matrix& M2) 
{
    return Matrix::add(M1,M2);
}

Matrix operator-(Matrix& M1, Matrix& M2) 
{
    return Matrix::subtract(M1,M2);
}

Matrix operator*(Matrix& M1, Matrix& M2) 
{
    return Matrix::multiply(M1,M2);
}

Matrix operator*(double  a, Matrix& M1) 
{   
    vector<vector<double>> temp;
    temp.resize(M1.getRows(), vector<double>(M1.getCols())); // Resize the values vector

    for(int i=0; i<M1.getRows(); i++)
    {
        for(int j=0; j<M1.getCols(); j++)
        {
            temp[i][j] = a*M1[i][j];
        }
    }
    Matrix M = temp;
    return M;
}

Matrix operator/(Matrix& M1, Matrix& M2) 
{   
    return Matrix::multiply(M1,M2.inverse());
}

Matrix operator^(Matrix& M1, int n) 
{   
    if(M1.getRows() != M1.getCols())
    {
        throw invalid_argument("Matrix must be square.");
    }

    Matrix M = Matrix::identity(M1.getRows());
    for(int i=1;i<=n; i++)
    {
        M = Matrix::multiply(M1,M);
    }
    return M;
}

// Output operator
std::ostream& operator<<(std::ostream& os, Matrix& M) 
{
    for (int i = 0; i < M.getRows(); ++i) {
        for (int j = 0; j < M.getCols(); ++j) {
            os << std::setw(10) << M[i][j] << " ";
        }
        os << std::endl;
    }
    return os;
}