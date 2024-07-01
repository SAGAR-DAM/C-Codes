#include<iostream>
#include<string>
#include<limits>
#include<vector>
#include<algorithm>
#include "D:\C++ codes\My_C++_modules\np.hpp"

using namespace std;

int main()
{

    // Example usage
    complex a(1.0, 2.0);
    complex b(3.0, 4.0);

    // std::cout << "complex a: " << a << std::endl;
    // std::cout << "complex b: " << b << std::endl;

    // Exponentiation
    complex c = a ^ b;
    std::cout << "a ^ b = " <<c;

    // complex c = complex(2,-5);
    // complex c1 = complex::acsccomplex(c);
    // cout<<"acsc(c): ";
    // c1.getcomplex();
    // c1.getcomplex_polar();


    // int x = printf("hello");
    // cout<<x;

    // string y = "dkajdkajd";
    // cout<<endl<<y.length()<<endl;  

    vector<double> array;
    for (int i=1; i<=100; i++)
    {
        array.push_back(i);
    }
    // for(int i=0; i<array.size(); i++)
    // {
    //     cout<<array[i]<<" ";
    // }
    // cout<<endl;

    // cout<<"max (array[1:6]): "<<*(std::max_element(array.begin()+1,array.begin()+6))<<endl;
    // cout<<"index of value 25.2: "<<(std::find(array.begin(),array.end())<<endl,25.2);

    // cout<<find_index(array,25.2)<<endl;

    // vector<double> newarray = cutarr(array,20,50);
    // for(int i=0; i<newarray.size(); i++) { cout<<newarray[i]<<" ";}


    // cout<<"mean array[20:50]: "<<meanarr(cutarr(array,20,50));
    // cout<<"stdev array[20:50]: "<<stdarr(cutarr(array,20,50));

    // vector<double> newarr = cutarr(array, find_closest_index(array,20.1),find_closest_index(array,45.2));
    // cout<<"newarr[0]: "<<newarr[0]<<" newarr[-1]: "<<newarr[newarr.size()-1]<<endl;


    std::cin.get();
    return(0);
}