// np.h

#ifndef NP_H  // Include guards to prevent multiple inclusion
#define NP_H

#include<iostream>
#include<string>
#include<limits>
#include<vector>
#include<algorithm>
#include <cmath>

using namespace std;

// OPERATOR DECLARATION



// FUNCTION DECLARATION:
int find_index(const std::vector<double>& array, double value);
int find_closest_index(const std::vector<double>& array, double value);
std::vector<double> cutarr(std::vector<double> array, int minw, int maxw);
double meanarr(const std::vector<double>& vec);
double stdarr(const std::vector<double>& vec);







// CLASS DECLARATION:
class complex
{
private:
    double real, imag;
    double r = sqrt(real * real + imag * imag);
    float theta = atan2(imag, real);
    //const complex iota;

public:
    // Constructor
    complex(double r, double i);

    //void rect_to_polar();
    //void polar_to_rect();
    double getreal();
    double getimag();
    double get_r();
    float get_theta();
    void getcomplex();
    void getcomplex_polar();

    static double abscomplex(complex c);
    static float argcomplex(complex c);
    static complex addcomplex(complex c1, complex c2);
    static complex subtractcomplex(complex c1, complex c2);
    static complex multcomplex(complex c1, complex c2);
    static complex divcomplex(complex c1, complex c2);
    static complex complex_exponenciation(complex c, float n);
    static complex complexpow(complex c1, complex c2);
    static complex real_pow_complex(complex c, float n);
    static complex sincomplex(complex c);
    static complex coscomplex(complex c);
    static complex tancomplex(complex c);
    static complex csccomplex(complex c);
    static complex seccomplex(complex c);
    static complex cotcomplex(complex c);
    static complex logcomplex(complex c);
    static complex asincomplex(complex c);
    static complex acoscomplex(complex c);
    static complex atancomplex(complex c);
    static complex acsccomplex(complex c);
    static complex aseccomplex(complex c);
    static complex acotcomplex(complex c);

};

// Complex operator definition outside the class
complex operator+(complex& c1, complex& c2);
complex operator-(complex& c1, complex& c2);
complex operator*(complex& c1, complex& c2);
complex operator/(complex& c1, complex& c2);
complex operator^(complex& c1, complex& c2);
std::ostream& operator<<(std::ostream& os, complex& c);





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
        static bool isDiagonal(Matrix mat);
        static Matrix invertMatrix(Matrix A);
        static pair<Matrix, Matrix> qrDecomposition(Matrix A);
        static vector<double> eigenvalues(Matrix A, int maxIterations, double tolerance);
        static Matrix qrAlgorithm(Matrix A, int maxIterations, double tolerance);
        static Matrix similarity_transform(Matrix A);

};

Matrix operator+(Matrix& M1, Matrix& M2);
Matrix operator-(Matrix& M1, Matrix& M2);
Matrix operator*(Matrix& M1, Matrix& M2);
Matrix operator*(double a, Matrix& M);
Matrix operator/(Matrix& M1, Matrix& M2);
Matrix operator^(Matrix& M, int n);
std::ostream& operator<<(std::ostream& os,Matrix& M);

#endif // np.h
