#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <numeric>
#include <cmath>
#include <ctime>
#include <cstdio>  // for remove()

int main() {
    int samples, sum_length;

    // Take input from the user
    std::cout << "Enter the number of samples: ";
    std::cin >> samples;
    std::cout << "Enter the sample size (number of random numbers per sample): ";
    std::cin >> sum_length;

    std::vector<double> means;

    std::mt19937 gen(static_cast<unsigned>(std::time(nullptr)));
    std::uniform_real_distribution<> dis(0.0, 7.0); // x in [0, 7]

    // Generate sample means
    for (int i = 0; i < samples; ++i) {
        double sum = 0.0;
        for (int j = 0; j < sum_length; ++j) {
            double x = dis(gen);
            sum += 1.0 / (1.0 + x);  // Apply the transformation
        }
        means.push_back(sum / sum_length);  // Mean of transformed values
    }

    // Compute mean and std deviation
    double mean = std::accumulate(means.begin(), means.end(), 0.0) / means.size();
    double sq_sum = 0.0;
    for (double x : means) sq_sum += (x - mean) * (x - mean);
    double stddev = std::sqrt(sq_sum / means.size());

    // Standardize
    for (double& x : means) x = (x - mean) / stddev;

    // Write histogram data
    const char* filename = "data_hist.txt";
    std::ofstream outfile(filename);
    for (double x : means) {
        outfile << x << "\n";
    }
    outfile.close();

    // Write standard normal distribution points
    const char* normal_file = "normal_curve.txt";
    std::ofstream normal_out(normal_file);
    int num_points = 1000;
    for (int i = 0; i <= num_points; ++i) {
        double x = -4.0 + 8.0 * i / static_cast<double>(num_points); // linspace from -4 to 4
        double y = (1.0 / std::sqrt(2 * M_PI)) * std::exp(-0.5 * x * x);
        normal_out << x << " " << y << "\n";
    }
    normal_out.close();

    // Open GNUplot pipe
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title 'CLT: Standardized Sample Means from f(x)=1/(1+x)'\n");
        fprintf(gnuplotPipe, "set xlabel 'Value'\n");
        fprintf(gnuplotPipe, "set ylabel 'Density'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "set xrange [-4:4]\n");
        fprintf(gnuplotPipe, "bin_width = 0.1\n");
        fprintf(gnuplotPipe, "bin(x,width) = width*floor(x/width) + width/2.0\n");

        fprintf(gnuplotPipe,
            "plot '%s' using (bin($1,bin_width)):(1.0 / (%d * bin_width)) smooth freq "
            "with boxes lc rgb 'green' fill solid 1.0 border lc rgb 'green' title 'Sample Means (Green)', \\\n",
            filename, samples);

        fprintf(gnuplotPipe,
            "     '%s' using (bin($1,bin_width)):(1.0 / (%d * bin_width)) smooth freq "
            "with boxes lc rgb 'blue' fill transparent border lc rgb 'blue' notitle, \\\n",
            filename, samples);

        fprintf(gnuplotPipe,
            "     '%s' using 1:2 with lines lw 2 lc rgb 'red' title 'Standard Normal PDF'\n",
            normal_file);

        fflush(gnuplotPipe);

        std::cout << "Press Enter to close plot and delete files..." << std::endl;
        std::cin.ignore();  // Ignore newline from earlier input
        std::cin.get();     // Wait for user input

        // Close pipe
        pclose(gnuplotPipe);

        // Clean up files
        if (remove(filename) != 0) {
            std::cerr << "Error deleting " << filename << std::endl;
        } else {
            std::cout << "Deleted " << filename << std::endl;
        }

        if (remove(normal_file) != 0) {
            std::cerr << "Error deleting " << normal_file << std::endl;
        } else {
            std::cout << "Deleted " << normal_file << std::endl;
        }

    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }

    return 0;
}
