// np.h

#ifndef NP_H  // Include guards to prevent multiple inclusion
#define NP_H

#include<iostream>
#include<string>
#include<limits>
#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;

// OPERATOR DECLARATION



// FUNCTION DECLARATION:
double factorial(int n);
int find_index(const std::vector<double>& array, double value);
int find_closest_index(const std::vector<double>& array, double value);
std::vector<double> cutarr(std::vector<double> array, int minw, int maxw);
double meanarr(const std::vector<double>& vec);
double stdarr(const std::vector<double>& vec);
std::vector<double> generateRandomdoubleVector(int size, double min_value, double max_value);
std::vector<int> generateRandomintVector(int size, int min_value, int max_value) ;
// Matrix generateRandomIntMatrix(int rows, int cols, int min_value, int max_value)  ;
// std::vector<std::vector<double>> generateRandomDoubleMatrix(int rows, int cols, double min_value, double max_value) ;






// CLASS DECLARATION:
class Complex
{
private:
    double real, imag;
    double r = sqrt(real * real + imag * imag);
    float theta = atan2(imag, real);
    //const Complex iota;

public:
    // Constructor
    Complex(double r, double i);

    //void rect_to_polar();
    //void polar_to_rect();
    double getreal();
    double getimag();
    double get_r();
    float get_theta();
    void getcomplex();
    void getcomplex_polar();
    Complex bar();

    static Complex conjugate(Complex c);
    static double abscomplex(Complex c);
    static float argcomplex(Complex c);
    static Complex addcomplex(Complex c1, Complex c2);
    static Complex subtractcomplex(Complex c1, Complex c2);
    static Complex multcomplex(Complex c1, Complex c2);
    static Complex divcomplex(Complex c1, Complex c2);
    static Complex complex_exponenciation(Complex c, float n);
    static Complex complexpow(Complex c1, Complex c2);
    static Complex real_pow_complex(Complex c, float n);
    static Complex sincomplex(Complex c);
    static Complex coscomplex(Complex c);
    static Complex tancomplex(Complex c);
    static Complex csccomplex(Complex c);
    static Complex seccomplex(Complex c);
    static Complex cotcomplex(Complex c);
    static Complex logcomplex(Complex c);
    static Complex asincomplex(Complex c);
    static Complex acoscomplex(Complex c);
    static Complex atancomplex(Complex c);
    static Complex acsccomplex(Complex c);
    static Complex aseccomplex(Complex c);
    static Complex acotcomplex(Complex c);

};

// Complex operator definition outside the class
Complex operator+(Complex c1, Complex c2);
Complex operator-(Complex c1, Complex c2);
Complex operator*(Complex c1, Complex c2);
Complex operator/(Complex c1, Complex c2);
Complex operator^(Complex c1, Complex c2);
std::ostream& operator<<(std::ostream& os, Complex& c);





class Matrix 
{
    private:
        vector<vector<double>> data;
        int rows, cols;

    public:
        Matrix(int rows, int cols);
        Matrix(vector<vector<double>>& values);

        int getRows();
        int getCols();

        vector<double>& operator[](int i);
        //vector<double>& operator[](int i);
        void print();

        Matrix transpose();
        double determinant();
        Matrix inverse();

        static Matrix add(Matrix a, Matrix b);
        static Matrix subtract(Matrix A, Matrix B);
        static Matrix multiply(Matrix a, Matrix b);
        static Matrix identity(int size);
        static double trace(Matrix M);
        static bool isDiagonal(Matrix mat);
        static Matrix invertMatrix(Matrix A);
        static pair<Matrix, Matrix> qrDecomposition(Matrix A);
        static vector<double> eigenvalues(Matrix A, int maxIterations, double tolerance);
        static Matrix qrAlgorithm(Matrix A, int maxIterations, double tolerance);
        static Matrix similarity_transform(Matrix A);
        static Matrix exp_matrix(Matrix M);
        static Matrix cos_matrix(Matrix M);
        static Matrix sin_matrix(Matrix M);
        static Matrix special_matrix(Matrix M, int k);
        static vector<double> coeff_charac(Matrix M);
};


Matrix operator+(Matrix M1, Matrix M2);
Matrix operator+(Matrix M,double  a);
Matrix operator-(Matrix M1, Matrix M2);
Matrix operator-(Matrix M,double  a);
Matrix operator*(Matrix M1, Matrix M2);
Matrix operator*(double a, Matrix M);
Matrix operator/(Matrix M1, Matrix M2);
Matrix operator/(double a, Matrix M);
Matrix operator^(Matrix M, int n);
std::ostream& operator<<(std::ostream& os,Matrix& M);
Matrix convertIntVector2ToMatrix(const std::vector<std::vector<int>>& int_matrix);
Matrix getRandomIntMatrix(int rows, int cols, int min_value, int max_value) ;
Matrix convertDoubleVector2ToMatrix(std::vector<std::vector<double>>& double_matrix);
Matrix getRandomDoubleMatrix(int rows, int cols, double min_value, double max_value);
Matrix convertIntToDoubleMatrix(Matrix& int_matrix) ;

#endif // np.h
