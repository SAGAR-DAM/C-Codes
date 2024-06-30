#include <iostream>
#include "D:\C++ codes\My_C++_modules\complex.hpp"
#include <fstream>
#include <cmath>
#include <cstdio> // For remove function
using namespace std;


complex repeatpower(complex c, int n)
{
    if(n<0)
    {
        cout<<"Wrong input; n should not be negative!!"<<endl;
    }
    if(n==0) 
    {
        complex c1 = complex(1,0);
        return c1;
    }
    else 
    {
        complex c1 = c.complexpow(c,repeatpower(c,n-1));
        return c1;
    }
};



int main() {
    complex z=complex(0,1.4446678);
    complex temp_complex = complex(1,0);

    // Generate data for y = sin(x)
    std::ofstream data1("data_powertower.dat");
    for (int n = 0; n <= 150; n += 1) {
        temp_complex = repeatpower(z,n);
        double x = temp_complex.getreal();
        double y = temp_complex.getimag();
        data1 << x << " " << y << "\n";
    }
    data1.close();

    // Open a pipe to GNUPLOT
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    //cout<<bool(gnuplotPipe)<<endl;
    if (gnuplotPipe) {
        // Send plotting commands to GNUPLOT
        fprintf(gnuplotPipe, "set title 'Powertower (z^{z^{z^{...^{z}}}}|_{n times}) plot of z = %.4f + %.4fi'\n", z.getreal(), z.getimag());
        fprintf(gnuplotPipe, "set xlabel 'Re(z)'\n");
        fprintf(gnuplotPipe, "set ylabel 'Im(z)'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "plot 'data_powertower.dat' title 'power towers' with linespoints pointtype 7 pointsize 1 linecolor rgb 'red' linewidth 1, ");
        fprintf(gnuplotPipe, "'data_powertower.dat' with points pointtype 7 pointsize 2 linecolor rgb 'blue'\n");
        fflush(gnuplotPipe);

        // Wait for the user to close the plot window
        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();

        // Close the pipe to GNUPLOT
        pclose(gnuplotPipe);

        // Delete the data files
        if (remove("data_powertower.dat") != 0) {
            std::cerr << "Error deleting data_powertower.dat" << std::endl;
        } else {
            std::cout << "Deleted data_powertower.dat" << std::endl;
        }
    } 
    
    else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }

    return 0;
}