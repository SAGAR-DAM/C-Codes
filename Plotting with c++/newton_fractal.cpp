#include <iostream>
#include <complex>
#include <vector>
#include <random>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <limits>

using namespace std;
using namespace chrono;

const int ITERATION = 50;
const int RES = 200;
const int rootnumber = 4;

// Function to create random roots
vector<complex<double>> create_random_root(int n) 
{
    vector<complex<double>> root;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(-1.0, 1.0);

    for (int i = 0; i < n; ++i) 
    {
        complex<double> z(dis(gen), dis(gen));
        root.push_back(z);
    }
    return root;
}

// Function to create n-th roots of unity
vector<complex<double>> create_circ_root(int n) 
{
    vector<complex<double>> roots;
    double angle_increment = 2.0 * M_PI / n;

    for (int k = 0; k < n; ++k) 
    {
        double angle = k * angle_increment;
        complex<double> z(cos(angle), sin(angle));
        roots.push_back(z);
    }

    return roots;
}

// Function to compute f(z)
complex<double> f(complex<double>& z, const vector<complex<double>>& root) {
    complex<double> val = 1.0;
    for (const auto& r : root) {
        val *= (z - r);
    }
    return val;
}

// Function to compute df(z)
complex<double> df(complex<double>& z, const vector<complex<double>>& root) {
    complex<double> val = 0.0;
    for (size_t i = 0; i < root.size(); ++i) 
    {
        complex<double> mult = 1.0;
        for (size_t j = 0; j < root.size(); ++j) 
        {
            if (j != i) {
                mult *= (z - root[j]);
            }
        }
        val += mult;
    }
    return val;
}

// Function to perform a single iteration step
complex<double> iteration_step(complex<double>& z, const vector<complex<double>>& root) 
{
    if(abs(f(z,root))<=0.000001 || abs(df(z,root))<=0.000001)
    {
        return z;
    }
    return z - f(z, root) / df(z, root);
}

int main() {
    auto t_start = high_resolution_clock::now();

    // Generate random roots
    //vector<complex<double>> root = create_random_root(rootnumber);
    vector<complex<double>> root = create_circ_root(rootnumber);

    cout << "Roots:" << endl;
    for (size_t i = 0; i < root.size(); i++) 
    {
        cout << "z" << i + 1 << ": " << root[i] << endl;
    }

    // Define the range for the plot
    double maxReal = 1.1 * max(abs(real(root[0])), abs(imag(root[0])));
    for (const auto& r : root) 
    {
        maxReal = max(maxReal, 1.1 * max(abs(real(r)), abs(imag(r))));
    }

    vector<double> x(2 * RES * (maxReal + 1));
    vector<double> y(2 * RES * (maxReal + 1));

    for (int i = 0; i < x.size(); ++i) {
        x[i] = -maxReal + (2 * maxReal / (x.size() - 1)) * i;
    }
    for (int i = 0; i < y.size(); ++i) {
        y[i] = -maxReal + (2 * maxReal / (y.size() - 1)) * i;
    }

    vector<complex<double>> z(x.size() * y.size());
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < y.size(); j++) {
            z[i * x.size() + j] = complex<double>(x[i], y[j]);
        }
    }

    // for(int i=0; i<z.size(); i++)
    // {
    //     cout<<real(z[i])<<" "<<imag(z[i])<<"\n";
    // }

    for (int i = 0; i < ITERATION; i++) 
    {
        for (size_t j = 0; j < z.size(); j++) 
        {
            z[j] = iteration_step(z[j], root);
        }
    }

    // Find min and max values of |f(z)| for proper color scaling
    double minVal = numeric_limits<double>::max();
    double maxVal = numeric_limits<double>::lowest();
    vector<int> f_values(z.size());

    int minindex;
    double dist, val;
    for (size_t i = 0; i < z.size(); i++) 
    {
        minindex = 0;
        dist = abs(z[i]-root[0]);
        for(int j=1; j<root.size(); j++)
        {
            val = abs(z[i] - root[j]);
            if(val<=dist)
            {
                minindex = j;
                dist = val;
            }
        }
        f_values[i] = minindex;
        val = abs(f(z[i], root));
        if (val < minVal) minVal = val;
        if (val > maxVal) maxVal = val;
    }

    // Avoid empty color bar range
    if (minVal == maxVal) {
        minVal -= 1.0;
        maxVal += 1.0;
    }

    // Create data file for Gnuplot
    ofstream dataFile("data.txt");
    // for (int i = 0; i < y.size(); ++i) {
    //     for (int j = 0; j < x.size(); ++j) {
    //         dataFile << real(z[i * x.size() + j]) << " " << imag(z[i * x.size() + j]) << " " << f_values[i * x.size() + j] << "\n";
    //     }
    //     dataFile << "\n";
    // }
    for (int i = 0; i < x.size(); i++)
    {
        for(int j=0; j < y.size(); j++)
        {
            dataFile << x[i] << " " << y[j] << " " << f_values[i * x.size() + j] << "\n";      
        } 
    }
    dataFile << "\n";
    dataFile.close();

    // Create Gnuplot script
    std::ofstream gnuplotScript("plot.gnuplot");
    // gnuplotScript << "set term wxt\n";  // Change to 'qt' if you prefer Qt terminal
    gnuplotScript << "set term wxt size 900,900\n";  // Adjust size as needed
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

    gnuplotScript << "set xlabel 'Re(z)'\n";
    gnuplotScript << "set ylabel 'Im(z)'\n";
    gnuplotScript << "set cblabel 'Root Index'\n";
    gnuplotScript << "set xtics font ',10'\n"; // Set smaller font for ticks and labels
    gnuplotScript << "set ytics font ',10'\n";
    gnuplotScript << "set xlabel 'Re(z)' font ',10'\n";
    gnuplotScript << "set ylabel 'Im(z)' font ',10'\n";
    gnuplotScript << "set cblabel 'Root Index' font ',10'\n";
    gnuplotScript << "set cbtics font ',10'\n";
    gnuplotScript << "set xrange [" << -maxReal << ":" << maxReal << "]\n";
    gnuplotScript << "set yrange [" << -maxReal << ":" << maxReal << "]\n";
    gnuplotScript << "set title 'Newton Fractal' font ',10'\n";
    gnuplotScript << "unset key\n"; // Disable the plot key (legend)
    gnuplotScript << "splot 'data.txt' using 1:2:3 with image\n";
    gnuplotScript.close();

    // Execute Gnuplot script
    system("gnuplot -persist plot.gnuplot");

    cout << "Interactive colormap displayed." << endl;
    // Wait for the user to close the plot window
    cout << "Press enter to exit..." << endl;
    cin.get();

    auto t_end = high_resolution_clock::now();
    duration<double> elapsed = t_end - t_start;
    cout << "\n\nImage shape: " << x.size() << "x" << y.size() << endl;
    cout << "Time taken: " << elapsed.count() << " seconds" << endl;

    // Delete the data files
    if (remove("data.txt") != 0) 
    {
        std::cerr << "Error deleting data.txt" << std::endl;
    } 
    else 
    {
        std::cout << "Deleted data.txt" << std::endl;
    }

    return 0;
}
