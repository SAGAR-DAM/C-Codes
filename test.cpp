#include<iostream>
#include<cmath>
#include<limits>
using namespace std;

int main()
{
    int n;
    int sum=0;
    int i,j;

    cout<<"Enter n:  ";
    cin>>n;
    cout<<"\n";

    for(int i=0; i<=n; i++)
    {
        sum = sum+i;
    }

    for(int i=0; i<=n; i++)
    {
        for(int j=1; j<=i; j++)
        {
            cout<<"*";
        }
        cout<<"  "<<i<<endl;
    }

    cout<<"The sum: "<<sum<<endl;
    cout<<"Direct sum:  "<<n*(n+1)/2;
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cin.get();
    return 0;
}