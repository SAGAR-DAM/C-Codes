#include "D:\C++ codes\My_C++_modules\complex.hpp"

using namespace std;

// Constructor
complex::complex(double r, double i) : real(r), imag(i) {}

void complex::rect_to_polar()
{
    r = sqrt(real * real + imag * imag);
    theta = atan2(imag, real);
}

void complex::polar_to_rect()
{
    real = r * cos(theta);
    imag = r * sin(theta);
}

double complex::getreal()
{
    return real;
}

double complex::getimag()
{
    return imag;
}

double complex::get_r(complex c)
{
    c.rect_to_polar();
    return r;
}

float complex::get_theta(complex c)
{
    c.rect_to_polar();
    return theta;
}

void complex::getcomplex()
{
    if (imag >= 0)
    {
        cout << "The Complex number is: " << real << " + " << imag << "j" << endl;
    }
    else
    {
        cout << "The Complex number is: " << real << " " << imag << "j" << endl;
    }
}

void complex::getcomplex_polar()
{
    rect_to_polar();
    cout << "The complex number in polar form: " << r << " exp(" << theta << "j)" << endl;
}

double complex::abscomplex(complex c)
{
    return pow((c.getreal()*c.getreal()+c.getimag()*c.getimag()),0.5);
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
    c1.rect_to_polar();
    c2.rect_to_polar();
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
    c1.rect_to_polar();
    c2.rect_to_polar();
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
    c.rect_to_polar();
    double r_exp = pow(c.r, n);
    float theta_exp = n * c.theta;

    double real_exp = r_exp * cos(theta_exp);
    double imag_exp = r_exp * sin(theta_exp);
    return complex(real_exp, imag_exp);
}

complex complex::complexpow(complex c1, complex c2)
{
    c1.rect_to_polar();
    c2.rect_to_polar();
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
    c_sin.rect_to_polar();
    return c_sin;
}

complex complex::coscomplex(complex c)
{
    double real_cos = cos(c.real)*cosh(c.imag), imag_cos = -sin(c.real)*sinh(c.imag);
    complex c_cos = complex(real_cos, imag_cos);
    c_cos.rect_to_polar();
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
    c.rect_to_polar();
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