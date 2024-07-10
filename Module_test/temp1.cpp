#include<iostream>
#include<string>
#include<limits>
#include"D:\C++ codes\My_C++_modules\np.hpp"
#include<cmath>
#include<algorithm>
#include<complex>
//#include<np.hpp>

using namespace std;
//using Matrix::Matrix;

// double factorial(int n)
// {
//     if(n==1)
//     {
//         double a = 1;
//         return a;
//     }

//     else
//     {
//         double val =  n*factorial(n-1);
//         return val;
//     }
// }

// Matrix exp_matrix(Matrix M)
// {
//     Matrix M1 = Matrix::identity(M.getRows());
//     for(int i=1;i<=100;i++)
//     {
//         M1 = M1 + 1/(factorial(i))*(M^i);
//     }
//     return M1;
// }

// Matrix cos_matrix(Matrix M)
// {
//     Matrix M1 = Matrix::identity(M.getRows());
//     for(int i=1; i<=100; i++)
//     {
//         M1 = M1 + pow(-1,i)/factorial(2*i)*(M^(2*i));
//     }
//     return M1;
// }

// Matrix sin_matrix(Matrix M)
// {
//     Matrix M1 = M;
//     for(int i=1; i<=100; i++)
//     {
//         M1 = M1 + pow(-1,i)/factorial(2*i+1)*(M^(2*i+1));
//     }
//     return M1;
// }

// // vector<double> coeff_charac(Matrix M)
// // {
// //     if (M.getRows() != M.getCols()) 
// //     {
// //         throw invalid_argument("Matrix must be square !!");
// //     }

// //     int n = M.getRows();
// //     vector<double> coeff;

// //     if(n==1)
// //     {
// //         coeff.push_back(-Matrix::trace(M));

// //     }
// //     else if(n==2)
// //     {
// //         coeff.push_back(M.determinant());
// //         coeff.push_back(-Matrix::trace(M));
// //     }
// //     else if(n==3)
// //     {
// //         coeff.push_back(M.determinant());
// //         coeff.push_back(1/2*(pow(Matrix::trace(M),2)-Matrix::trace(M^2)));
// //         coeff.push_back(-Matrix::trace(M));
// //     }
// //     else if(n==4)
// //     {
// //         coeff.push_back(M.determinant());
// //         coeff.push_back(-1/6*(pow(Matrix::trace(M),3)-2*Matrix::trace(M)*Matrix::trace(M^2)+2*Matrix::trace(M^3)));
// //         coeff.push_back(1/2*(pow(Matrix::trace(M),2)-Matrix::trace(M^2)));
// //         coeff.push_back(-Matrix::trace(M));
// //     }
// //     // else
// //     // {
// //     //     double ck;
// //     //     coeff.push_back(pow(-1,n)*M.determinant());
// //     //     for(int k=1;k<=n-1;k++)
// //     //     {
// //     //         ck= -1/k*Matrix::trace(M^k);
// //     //         for(int i=1; i<=k-1; i++)
// //     //         {
// //     //             ck = ck - 1/k*Matrix::trace(M^i)*coeff[k-i];
// //     //             cout<<i<<endl;
// //     //         }
// //     //         // ck = ck*pow(-1,k+1)/k;
// //     //         coeff.push_back(ck);
// //     //     }
// //     // }

// //     coeff.push_back(1);
// //     // Reverse the vector
// //     std::reverse(coeff.begin(), coeff.end());
// //     return coeff;
// // }


// Matrix special_matrix(Matrix M, int k)
// {
//     Matrix A = Matrix::identity(k);
//     for(int i=0; i<k; i++)
//     {
//         for(int j=0; j<k; j++)
//         {
//             if(j<=i)
//             {
//                 A[i][j]=Matrix::trace(M^(i-j+1));
//             }
//             if(j==i+1)
//             {
//                 A[i][j]=k-j;
//             }
//         }
//     }
//     return A;
// }

// vector<double> coeff_charac(Matrix M)
// {
//     if (M.getRows() != M.getCols()) 
//     {
//         throw invalid_argument("Matrix must be square !!");
//     }
//     int n = M.getRows();
//     vector<double> coeff;
//     for(int i = 0; i<=n; i++)
//     {
//         if(i==0)
//         {
//             coeff.push_back(1);
//         }
//         else
//         {
//             Matrix A = special_matrix(M,i);
//             double val = pow(-1,i)*A.determinant()/factorial(i);
//             coeff.push_back(val);
//         }
//     }
//     //std::reverse(coeff.begin(), coeff.end());
//     return coeff;
// }

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
    
    // = {
    //     {4, 7, 2},
    //     {3, 6, 1},
    //     {2, 5, 1}
    // };

    // = {
    //     {3, 0},
    //     {0, 1}
    // };

    // = {
    //     {4, 1, 0, 0},
    //     {2, 6, 0, 0},
    //     {0, 0, 1, 0},
    //     {0, 0, 0, 7}
    // };

    Matrix mat = values;
    cout<<"trace: "<<Matrix::trace(mat)<<endl;
    cout<<"det: "<<mat.determinant()<<endl;
    cout<<"term: "<<Matrix::trace(mat^2)<<endl;
    // cout<<mat;
    // double a = 5;
    // Matrix M = a*mat;
    // cout<<M;

    // M = mat*mat;
    // cout<<M;
    // Matrix M = sin_matrix(mat);
    // cout<<M<<endl;
    // double tr = Matrix::trace(M);
    // cout<<tr;

    vector<double> ck = Matrix::coeff_charac(mat);
    for(int i=0; i<ck.size(); i++)
    {
        cout<<ck[i]<<", ";
    }

    // Matrix A = special_matrix(mat,0);
    // cout<<A;
    std::cin.get();
    return(0);
}