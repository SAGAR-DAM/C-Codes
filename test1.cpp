#include<iostream>
#include<limits>
using namespace std;

long double fact(int n)
{
    if(n == 1 || n == 0)
    {
        return 1;
    }
    // else
    // {
    //     return n * fact(n-1);
    // }
    return n*fact(n-1);
}

int main()
{
    int n;
    cout<<"Enter the number: ";
    cin>>n;

    if(n < 0)
    {
        cout << "Factorial is not defined for negative numbers." << endl;
    }
    else
    {
        long double result = fact(n);
        cout << n << "! : " << result;
    }
    cout<<endl;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cin.get();
    return 0;
}