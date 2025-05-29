#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio> // For popen function

int main() {
    // Generate data for 3D line plot (x, y, z)
    std::ofstream data("3d_data.dat");
    for (double t = 0; t <= 10; t += 0.1) {
        double x = sin(t);
        double y = cos(t);
        double z = t;  // This creates a spiral
        data << x << " " << y << " " << z << "\n";
    }
    data.close();

    // Write GNUPLOT script for 3D line plot
    std::ofstream gp_script("plot_3d.gp");
    gp_script << "set title '3D Line Plot'\n";
    gp_script << "set xlabel 'x'\n";
    gp_script << "set ylabel 'y'\n";
    gp_script << "set zlabel 'z'\n";
    gp_script << "set grid\n";
    gp_script << "set terminal wxt size 800,600\n";  // Interactive window
    gp_script << "splot '3d_data.dat' with lines title '3D Spiral'\n";  // Plotting the 3D line
    gp_script.close();

    // Execute GNUPLOT script
    FILE *gnuplotPipe = popen("gnuplot", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "load 'plot_3d.gp'\n");
        fflush(gnuplotPipe);

        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();

        pclose(gnuplotPipe);
    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }

    // Optionally, delete the data file
    if (remove("3d_data.dat") != 0) {
        std::cerr << "Error deleting 3d_data.dat" << std::endl;
    } else {
        std::cout << "Deleted 3d_data.dat" << std::endl;
    }

    if (remove("plot_3d.gp") != 0) {
        std::cerr << "Error deleting plot_3d.gp" << std::endl;
    } else {
        std::cout << "Deleted plot_3d.gp" << std::endl;
    }

    return 0;
}
