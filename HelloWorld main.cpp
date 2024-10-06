#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include "D:\C++ codes\My_C++_modules\np.hpp"
#include <limits>  
#include <regex>
#include <chrono>   // to create new random numbers each time. Visit: https://chatgpt.com/share/6700cacc-8640-8013-bc41-94054ab0dc0b

using namespace std;

int main()
{
    Matrix A=getRandomIntMatrix(30,30,-9,9);
    std::cout<<"Hello World!!\n";
    cout<<A<<endl;
    return 0;
}