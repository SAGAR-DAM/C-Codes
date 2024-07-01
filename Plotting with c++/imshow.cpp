#include <iostream>
#include <fstream>
#include <cmath>

// Function to compute z = sin(x^2 + y^2) * exp(-(x^2 + y^2))
double computeZ(double x, double y) {
    // double r2 = x * x + y * y;
    // return sin(r2) / r2;

    double z=sqrt(x*x+y*y)+3*cos(sqrt(x*x+y*y))+5;

    return z;
}

int main() {
    int X_gridSize = 501, Y_gridSize = 501; // Define the resolution of the grid
    double xMin = -20, xMax = 20;
    double yMin = -20, yMax = 20;

    // Create a file to store the data
    std::ofstream dataFile("data.txt");

    // Generate data
    double x_step = (xMax - xMin) / X_gridSize, y_step = (yMax - yMin) / Y_gridSize;
    for (int i = 0; i <= X_gridSize; ++i) {
        double x = xMin + i * x_step;
        for (int j = 0; j <= Y_gridSize; ++j) {
            double y = yMin + j * y_step;
            double z = computeZ(x, y);
            dataFile << x << " " << y << " " << z << "\n";
        }
        dataFile << "\n";
    }
    dataFile.close();

    // Create Gnuplot script
    std::ofstream gnuplotScript("plot.gnuplot");
    gnuplotScript << "set term wxt\n";  // Change to 'qt' if you prefer Qt terminal
    gnuplotScript << "set view map\n";
    gnuplotScript << "set size ratio -1\n";

    // Define the "jet" colormap
    gnuplotScript << "set palette defined ( \\\n";
    gnuplotScript << "0 0 0 0, \\\n"; // Black
    gnuplotScript << "0.1 0 0 0.5, \\\n"; // Dark Blue
    gnuplotScript << "0.2 0 0 1, \\\n"; // Blue
    gnuplotScript << "0.3 0 0.5 1, \\\n"; // Cyan Blue
    gnuplotScript << "0.4 0 1 1, \\\n"; // Cyan
    gnuplotScript << "0.5 0.5 1 0.5, \\\n"; // Light Green
    gnuplotScript << "0.6 1 1 0, \\\n"; // Yellow
    gnuplotScript << "0.7 1 0.5 0, \\\n"; // Orange
    gnuplotScript << "0.8 1 0 0, \\\n"; // Red
    gnuplotScript << "0.9 0.5 0 0, \\\n"; // Dark Red
    gnuplotScript << "1 0.5 0 0 )\n"; // Dark Red


    gnuplotScript << "set xlabel 'x'\n";
    gnuplotScript << "set ylabel 'y'\n";
    gnuplotScript << "set cblabel 'z'\n";
    gnuplotScript << "set xtics font ',10'\n"; // Set smaller font for ticks and labels
    gnuplotScript << "set ytics font ',10'\n";
    gnuplotScript << "set xlabel 'x' font ',10'\n";
    gnuplotScript << "set ylabel 'y' font ',10'\n";
    gnuplotScript << "set cblabel 'z' font ',10'\n";
    gnuplotScript << "set cbtics font ',10'\n";
    gnuplotScript << "set xrange [" << xMin << ":" << xMax << "]\n";
    gnuplotScript << "set yrange [" << yMin << ":" << yMax << "]\n";
    gnuplotScript << "set title 'contour plot' font ',10'\n";
    gnuplotScript << "unset key\n"; // Disable the plot key (legend)
    gnuplotScript << "splot 'data.txt' using 1:2:3 with image\n";
    gnuplotScript.close();

    // Execute Gnuplot script
    system("gnuplot -persist plot.gnuplot");

    std::cout << "Interactive colormap displayed." << std::endl;
    // Wait for the user to close the plot window
    std::cout << "Press enter to exit..." << std::endl;
    std::cin.get();

    // Delete the data files
    if (remove("data.txt") != 0) {
        std::cerr << "Error deleting data.txt" << std::endl;
    } else {
        std::cout << "Deleted data.txt" << std::endl;
    }
    return 0;
}
