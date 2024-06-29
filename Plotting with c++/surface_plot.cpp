#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio> // For popen function

int main() {
    // Generate surface data
    std::ofstream data("surface_data.dat");
    for (double x = -5.0; x <= 5.0; x += 0.1) {
        for (double y = -5.0; y <= 5.0; y += 0.1) {
            double z = exp(-x*x-y*y)*cos(0.2*(x*x+y*y));
            data << x << " " << y << " " << z << "\n";
        }
        data << "\n";
    }
    data.close();

    // Write GNUPLOT script
    std::ofstream gp_script("plot_surface.gp");
    gp_script << "set title 'Interactive Surface Plot'\n";
    gp_script << "set xlabel 'x'\n";
    gp_script << "set ylabel 'y'\n";
    gp_script << "set zlabel 'z'\n";
    gp_script << "set pm3d\n";
    gp_script << "set hidden3d\n";
    gp_script << "splot 'surface_data.dat' with lines\n";
    gp_script.close();

    // Execute GNUPLOT script
    FILE *gnuplotPipe = popen("gnuplot", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "load 'plot_surface.gp'\n");
        fflush(gnuplotPipe);

        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();

        pclose(gnuplotPipe);
    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }

    // Optionally, delete the data file
    if (remove("surface_data.dat") != 0) {
        std::cerr << "Error deleting surface_data.dat" << std::endl;
    } else {
        std::cout << "Deleted surface_data.dat" << std::endl;
    }

    if (remove("plot_surface.gp") != 0) {
        std::cerr << "Error deleting plot_surface.gp" << std::endl;
    } else {
        std::cout << "Deleted plot_surface.gp" << std::endl;
    }

    return 0;
}
