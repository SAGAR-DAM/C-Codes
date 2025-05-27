#include <iostream>
#include <string>
using namespace std;

struct Point
{
    int x=0, y=0;
    // Point(int x=0, int y=0)// default constructor
    // {
    //     this->x=x;
    //     this->y=y;
    // } 
    
};

int main()
{
    Point p{.y=3}, q{.x=5};
    cout<<"p.x: "<<p.x<<"; p.y: "<<p.y<<endl;
    cout<<"q.x: "<<q.x<<"; q.y: "<<q.y<<endl;
    return 0;
}