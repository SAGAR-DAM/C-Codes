#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio> // For remove function

int main() {
    // Generate data for y = sin(x)
    std::ofstream data1("data_sin.dat");
    for (double x = -20.0; x <= 20.0; x += 0.01) {
        double y = sin(x);
        data1 << x << " " << y << "\n";
    }
    data1.close();

    // Generate data for y = exp(cos(x))
    std::ofstream data2("data_exp_cos.dat");
    for (double x = -20.0; x <= 20.0; x += 0.01) {
        double y = exp(cos(x));
        data2 << x << " " << y << "\n";
    }
    data2.close();
    
    // Generate data for y = cos(x)
    std::ofstream data3("data_cos.dat");
    for (double x = -20.0; x <= 20.0; x += 0.01) {
        double y = cos(x);
        data3 << x << " " << y << "\n";
    }
    data3.close();

    // Open a pipe to GNUPLOT
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        // Send plotting commands to GNUPLOT
        fprintf(gnuplotPipe, "set title 'Multiplot Example'\n");
        fprintf(gnuplotPipe, "set xlabel 'x'\n");
        fprintf(gnuplotPipe, "set ylabel 'y'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "plot 'data_sin.dat' title 'y = sin(x)' with lines linecolor rgb 'green' linewidth 2, ");
        fprintf(gnuplotPipe, "'data_cos.dat' title 'y = cos(x)' with lines linecolor rgb 'black' linewidth 2, ");
        fprintf(gnuplotPipe, "'data_exp_cos.dat' title 'y = exp(cos(x))' with linespoints pointtype 7 pointsize 1 linecolor rgb 'red' linewidth 2\n");
        fflush(gnuplotPipe);

        // Wait for the user to close the plot window
        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();

        // Close the pipe to GNUPLOT
        pclose(gnuplotPipe);

        // Delete the data files
        if (remove("data_sin.dat") != 0) {
            std::cerr << "Error deleting data_sin.dat" << std::endl;
        } else {
            std::cout << "Deleted data_sin.dat" << std::endl;
        }

        if (remove("data_exp_cos.dat") != 0) {
            std::cerr << "Error deleting data_exp_cos.dat" << std::endl;
        } else {
            std::cout << "Deleted data_exp_cos.dat" << std::endl;
        }

        if (remove("data_cos.dat") != 0) {
            std::cerr << "Error deleting data_cos.dat" << std::endl;
        } else {
            std::cout << "Deleted data_cos.dat" << std::endl;
        }
    } 
    
    else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }

    return 0;
}
