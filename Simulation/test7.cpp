#include <iostream>
#include <vector>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>

using namespace std;

void plot_3d_scatter_segments(const std::vector<double>& posx, const std::vector<double>& posy, const std::vector<double>& posz, size_t n)
{
    // Step 1: Split data into segments of size n
    size_t total_points = posx.size();
    size_t num_segments = (total_points + n - 1) / n;

    // Write each segment to a separate file
    for (size_t s = 0; s < num_segments; ++s) {
        std::string filename = "segment_" + std::to_string(s) + ".dat";
        std::ofstream seg_file(filename);
        for (size_t i = s * n; i < std::min((s + 1) * n, total_points); ++i) {
            seg_file << posx[i] << " " << posy[i] << " " << posz[i] << "\n";
        }
        seg_file.close();
    }

    // Step 2: Create the Gnuplot script
    std::ofstream gp_script("plot_3d_scatter.gp");
    gp_script << "set title '3D Particle Trajectory (Scatter)'\n";
    gp_script << "set xlabel 'X Position (m)'\n";
    gp_script << "set ylabel 'Y Position (m)'\n";
    gp_script << "set zlabel 'Z Position (m)'\n";
    gp_script << "set grid\n";
    gp_script << "set terminal wxt size 800,600\n";
    gp_script << "splot ";

    // Add each segment to the splot command with different colors
    for (size_t s = 0; s < num_segments; ++s) {
        std::string filename = "segment_" + std::to_string(s) + ".dat";
        gp_script << "'" << filename << "' using 1:2:3 with points pointtype 7 pointsize 1.5 title 'Segment " << s << "'";
        if (s != num_segments - 1)
            gp_script << ", \\\n";
    }
    gp_script << "\n";
    gp_script.close();

    // Step 3: Execute the Gnuplot script
    FILE *gnuplotPipe = popen("gnuplot", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "load 'plot_3d_scatter.gp'\n");
        fflush(gnuplotPipe);

        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();

        pclose(gnuplotPipe);
    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }

    // Cleanup
    for (size_t s = 0; s < num_segments; ++s) {
        std::string filename = "segment_" + std::to_string(s) + ".dat";
        if (remove(filename.c_str()) != 0)
            std::cerr << "Error deleting " << filename << std::endl;
    }

    if (remove("plot_3d_scatter.gp") != 0)
        std::cerr << "Error deleting plot_3d_scatter.gp" << std::endl;
}


int main()
{
    // Generate sample 3D trajectory: a helix for visualization
    std::vector<double> posx;
    std::vector<double> posy;
    std::vector<double> posz;

    double t_max = 10.0;
    double dt = 0.1;
    for (double t = 0; t <= t_max; t += dt)
    {
        posx.push_back(cos(t));       // x = cos(t)
        posy.push_back(sin(t));       // y = sin(t)
        posz.push_back(t / (2 * M_PI));    // z = increasing linearly (1 full turn per 2π)
    }

    // Number of points per color segment
    size_t segment_length = 20;

    // Call plotting function
    plot_3d_scatter_segments(posx, posy, posz, segment_length);

    return 0;
}