#include<iostream>
#include<string>
#include<limits>
#include"D:\C++ codes\My_C++_modules\np.hpp"
#include<cmath>
#include<algorithm>
#include<complex>

using namespace std;


int main()
{
    vector<vector<double>> values
    = {
        {4, 7, 2, 5, 4, 6},
        {3, 6, 1, 6,-8,-6},
        {2, 5, 1,-5, 7, 1},
        {7,-7, 5, 9,-2,-1},
        {2, 5, 4, 1,-5, 7},
        {5, 1,-5, 6,-8,-6}
    };
    
    Matrix mat = values;
    vector<double> ck = Matrix::coeff_charac(mat);
    for (int i = 0; i < ck.size(); i++) {
        cout << ck[i] << " ";
    }
    cout << endl;
    return(0);
}