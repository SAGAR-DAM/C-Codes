#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <limits>
#include <cmath>
#include <algorithm>  // Include for std::min_element and std::distance

const int order = 3;
const int iteration = 20;
const int size = 201;
const double x_min = -1.0;
const double x_max = 1.0;
const double y_min = -1.0;
const double y_max = 1.0;

std::vector<std::complex<double>> roots;
std::vector<double> colors(order);

std::complex<double> f(std::complex<double> z) {
    return std::pow(z, order) - std::complex<double>(1.0, 0.0);
}

std::complex<double> df(std::complex<double> z) {
    return std::complex<double>(order, 0.0) * std::pow(z, order - 1);
}

void colorComplexPlane(std::vector<std::vector<double>>& output, std::vector<std::vector<std::complex<double>>>& z) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            std::vector<double> dist(order);
            for (int k = 0; k < order; ++k) {
                dist[k] = std::abs(z[i][j] - roots[k]);
            }
            auto minIter = std::min_element(dist.begin(), dist.end());
            int nearest = std::distance(dist.begin(), minIter);
            output[i][j] = colors[nearest];
        }
    }
}

int main() {
    // Generate roots
    for (int i = 0; i < order; ++i) {
        double angle = 2 * M_PI * i / order;
        roots.emplace_back(std::polar(1.0, angle));
    }

    // Generate colors
    for (int i = 0; i < order; ++i) {
        colors[i] = static_cast<double>(i) / order;
    }

    // Create grid
    std::vector<double> x(size), y(size);
    for (int i = 0; i < size; ++i) {
        x[i] = x_min + i * (x_max - x_min) / (size - 1);
        y[i] = y_min + i * (y_max - y_min) / (size - 1);
    }

    std::vector<std::vector<std::complex<double>>> z(size, std::vector<std::complex<double>>(size));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            z[i][j] = std::complex<double>(x[j], y[i]);
        }
    }

    std::vector<std::vector<double>> output(size, std::vector<double>(size, 0.0));

    // Iterate
    for (int it = 0; it < iteration; ++it) {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                z[i][j] = z[i][j] - f(z[i][j]) / df(z[i][j]);
            }
        }
    }

    // Color the complex plane
    colorComplexPlane(output, z);

    // Output data to file
    std::ofstream outFile("output.dat");
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            outFile << x[j] << " " << y[i] << " " << output[i][j] << "\n";
        }
        outFile << "\n";
    }
    outFile.close();

    // Create Gnuplot script
    std::ofstream gnuplotScript("plot.gnuplot");
    gnuplotScript << "set term wxt size 800, 800\n";  // Change to 'qt' if you prefer Qt terminal
    gnuplotScript << "set view map\n";
    gnuplotScript << "set size ratio -1\n";
    
    // Define the "jet" colormap
    gnuplotScript << "set palette defined ( \\\n";
    gnuplotScript << "0 0 0 0, \\\n";           // Black
    gnuplotScript << "0.1 0 0 0.5, \\\n";       // Dark Blue
    gnuplotScript << "0.2 0 0 1, \\\n";         // Blue
    gnuplotScript << "0.3 0 0.5 1, \\\n";       // Cyan Blue
    gnuplotScript << "0.4 0 1 1, \\\n";         // Cyan
    gnuplotScript << "0.5 0.5 1 0.5, \\\n";     // Light Green
    gnuplotScript << "0.6 1 1 0, \\\n";         // Yellow
    gnuplotScript << "0.7 1 0.5 0, \\\n";       // Orange
    gnuplotScript << "0.8 1 0 0, \\\n";         // Red
    gnuplotScript << "0.9 0.5 0 0, \\\n";       // Dark Red
    gnuplotScript << "1 0.5 0 0 )\n";           // Dark Red
    
    gnuplotScript << "set xlabel 'Re(z)'\n";
    gnuplotScript << "set ylabel 'Im(z)'\n";
    gnuplotScript << "set xtics font ',10'\n"; // Set smaller font for ticks and labels
    gnuplotScript << "set ytics font ',10'\n";
    gnuplotScript << "set xlabel 'Re(z)' font ',10'\n";
    gnuplotScript << "set ylabel 'Im(z)' font ',10'\n";
    gnuplotScript << "unset colorbox\n"; // Disable the colorbox
    gnuplotScript << "plot 'output.dat' using 1:2:3 with points pt 7 ps 0.5 palette notitle\n";
    gnuplotScript.close();

    return 0;
}
