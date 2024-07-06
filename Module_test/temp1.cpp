#include<iostream>
#include<string>
#include<limits>
#include"D:\C++ codes\My_C++_modules\np.hpp"

using namespace std;

int main()
{
    vector<vector<double>> values 
    = {
        {3,-2},
        {4,-1}
    };

    // = {
    //     {4, 1, 0, 0},
    //     {2, 6, 0, 0},
    //     {0, 0, 1, 0},
    //     {0, 0, 0, 7}
    // };

    // = {
    //     {4, 7, 2, 5},
    //     {3, 6, 1, 6},
    //     {2, 5, 1,-5},
    //     {7,-7, 5, 0}
    // };

    Matrix mat = values;
    cout<<mat;
    double a = 5;
    Matrix M = a*mat;
    cout<<M;

    // M = mat*mat;
    // cout<<M;

    std::cin.get();
    return(0);
}