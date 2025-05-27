#include <iostream>
#include <string>
using namespace std;

struct Complex
{
    int x, y;
    Complex(int x = 0, int y = 0) // default constructor
    {
        this->x = x;
        this->y = y;
    }

    struct __init__
    {
        int x = 0;
        int y = 0;
        Complex __build__() const { return Complex(x, y); }
    };
};

struct Real : public Complex {
    Real(int x = 0) : Complex(x, 0) {}  // only real part
};

struct Imag : public Complex {
    Imag(int y = 0) : Complex(0, y) {}  // only imaginary part
};

int main()
{
    Real a(3);
    a.y=6;
    cout<<"a: "<<a.x<<"+i"<<a.y<<endl;
    Imag b(2);
    cout<<"b: "<<b.x<<"+i"<<b.y<<endl;
    return 0;
}