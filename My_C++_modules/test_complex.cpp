#include <iostream>
#include "D:\C++ codes\My_C++_modules\complex.hpp"

int main()
{
    complex c1(3, 4);
    complex c2(1, 2);

    c1.getcomplex();
    c2.getcomplex();

    complex c3 = complex::addcomplex(c1, c2);
    c3.getcomplex();

    complex c4 = complex::multcomplex(c1, c2);
    c4.getcomplex();

    return 0;
}
