#include <iostream>
#include <string>
using namespace std;

struct Point
{
    int x, y;
    Point(int x = 0, int y = 0) // default constructor
    {
        this->x = x;
        this->y = y;
    }

    struct __init__
    {
        int x = 0;
        int y = 0;
        Point __build__() const { return Point(x, y); }
    };
};

struct function_param
{
    Point p;
    int n;
    int a;
    int b;

    function_param(Point p=Point(), int n = 0, int a = 0, int b = 0)
    {
        this->p = p;
        this->n = n;
        this->a = a;
        this->b = b;
    }

    struct __init__
    {
        Point p;
        int n = 0;
        int a = 0;
        int b = 0;
        function_param __build__() const { return function_param(p, n, a, b); }
    };
};


double function(function_param fp)
{
    double val = (fp.p.x + fp.p.y) * fp.n + fp.a + fp.b;
    return val;
}

double function(Point p=Point(), int n=0, int a=0, int b=0)
{
    double val = (p.x + p.y) * n + a + b;
    return val;
}


int main()
{
    Point p = Point::__init__{.y = 10}.__build__(); // x defaults to 0
    Point q(1, 2);
    cout << "p.x: " << p.x << "; p.y: " << p.y << endl;
    cout << "q.x: " << q.x << "; q.y: " << q.y << endl;

    cout<<"function val:  "<<function(function_param::__init__{.p=q,.n=1}.__build__())<<endl;
    cout<<"default function with struct val: "<<function(Point())<<endl;
    cout<<"function with parameters val: "<<function(q,2,0,1)<<endl;
    cout<<"default function with parameters val: "<<function();
    return 0;
}

