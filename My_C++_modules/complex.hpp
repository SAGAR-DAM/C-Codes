#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#include <iostream>
#include <cmath>

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

#endif // COMPLEX_HPP
