#include <iostream>
#include "D:\C++ codes\My_C++_modules\np.hpp"
using namespace std;

int main()
{
    // get a Complex (simplest program):
    //________________________________________
    // double a,b;
    // cout<<"real part: ";
    // cin>>a;

    // cout<<"real part: ";
    // cin>>b;

    // Complex c1=Complex(a,b);
    // c1.getcomplex();

    // Complex number in polar form:
    // Complex c=Complex(1,1);
    // c.getcomplex_polar();

    // Complex no abs and arg:
    // Complex c = Complex(1,1);
    // cout<<"abs: "<<c.get_r()<<" ; arg: "<<c.get_theta();

    // abs Complex:
    // Complex c=Complex(1,1);
    // double absc = Complex::abscomplex(c);
    // cout<<absc;

    // abs Complex:
    // Complex c=Complex(1,1);
    // double argc = Complex::argcomplex(c);
    // cout<<argc;

    // add two Complex number:
    //________________________________________
    // Complex c1 = Complex(2,4), c2 = Complex(5,-1);
    // Complex c3 = Complex::addcomplex(c1,c2);
    // c3.getcomplex();

    // subtract two Complex number:
    //________________________________________
    // Complex c1 = Complex(2,4), c2 = Complex(5,-1);
    // Complex c3 = Complex::subtractcomplex(c1,c2);
    // c3.getcomplex();

    // multiply two Complex number:
    // ________________________________________
    // Complex c1 = Complex(2,4), c2 = Complex(5,-1);
    // Complex c3 = Complex::multcomplex(c1,c2);
    // c3.getcomplex();
    // c3.getcomplex_polar();

    // division of two Complex number:
    //________________________________________
    // Complex c1 = Complex(2,4), c2 = Complex(5,-1);
    // Complex c3 = Complex::divcomplex(c1,c2);
    // c3.getcomplex();
    // c3.getcomplex_polar();

    // Complex to the power Complex number:
    // ________________________________________
    // Complex c1 = Complex(2,-5), c2 = Complex(2.6,-1);
    // Complex c3 = Complex::complexpow(c1,c2);
    // c3.getcomplex();
    // c3.getcomplex_polar();

    // real to the power Complex number:
    // ________________________________________
    // Complex c1 = Complex(2,-5);
    // float n = 2;
    // Complex c3 = Complex::real_pow_complex(c1,n);
    // c3.getcomplex();
    // c3.getcomplex_polar();

    // sin, cos, tan, ... of Complex number:
    // ________________________________________
    // Complex c = Complex(2,-5);
    // Complex c1 = Complex::sincomplex(c);
    // cout<<"sin(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    // c1 = Complex::coscomplex(c);
    // cout<<"cos(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    // c1 = Complex::tancomplex(c);
    // cout<<"tan(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    // c1 = Complex::csccomplex(c);
    // cout<<"csc(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    // c1 = Complex::seccomplex(c);
    // cout<<"sec(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    // c1 = Complex::cotcomplex(c);
    // cout<<"cot(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();


    // log of Complex number:
    // ________________________________________
    // Complex c = Complex(1.414,1.414);
    // Complex c1 = Complex::logcomplex(c);
    // cout<<"log(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    // asin, acos, atan, ... of Complex number:
    // ________________________________________
    // Complex c = Complex(2,-5);
    // Complex c1 = Complex::asincomplex(c);
    // cout<<"asin(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    // Complex c = Complex(2,-5);
    // Complex c1 = Complex::acoscomplex(c);
    // cout<<"acos(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    // Complex c = Complex(2,-5);
    // Complex c1 = Complex::atancomplex(c);
    // cout<<"atan(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    // Complex c = Complex(2,-5);
    // Complex c1 = Complex::acsccomplex(c);
    // cout<<"acsc(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    // Complex c = Complex(2,-5);
    // Complex c1 = Complex::aseccomplex(c);
    // cout<<"asec(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();

    Complex c = Complex(2,-5);
    Complex c1 = Complex::acotcomplex(c);
    cout<<"acot(c): ";
    c1.getcomplex();
    c1.getcomplex_polar();



    std::cin.get();
    return(0);
}