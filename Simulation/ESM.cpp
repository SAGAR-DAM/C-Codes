#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <limits>
// #include "D:\C++ codes\matplotlibcpp\include\matplotlibcpp.h"

// namespace plt = matplotlibcpp;


using namespace std;



double cm = 1e-2;
double mm = 1e-3;
double qe = 1.60217663e-19;
double me = 9.10938e-31;
double mH = 1.67262192e-27;
double kev_to_joule = 1.660217663e-16;
double c = 3e8;

double Ex = 0;
double Ey = 0;
double Ez = 0;
double Bx = -0.3;
double By = 0;
double Bz = 0;

double magnet_height = 20*mm;


struct Particle
{
    /* data */
    double x, y, z, vx, vy, vz, q, m, v, energy;
    Particle(double x = 0.0,
             double y = 0.0,
             double z = 0.0,
             double vx = 0.0,
             double vy = 0.0,
             double vz = 0.0,
             double q = 0.0,
             double m = 1.0)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->vx = vx;
        this->vy = vy;
        this->vz = vz;
        this->q = q;
        this->m = m;
        v = sqrt(vx * vx + vy * vy + vz * vz);
        energy = 0.5 * m * v * v;
    }

    struct __init__
    {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double vx = 0.0;
        double vy = 0.0;
        double vz = 0.0;
        double q = 0.0;
        double m = 1.0;

        Particle __build__() const
        {
            return Particle(x, y, z, vx, vy, vz, q, m);
        }
    };

    std::vector<double> posx;
    std::vector<double> posy;
    std::vector<double> posz;
};


void propagator(Particle &p, double t_max = 0, double dt = 0.001, int save_step=50)
{
    int steps = static_cast<int>(t_max / dt);
    for(int i=0; i<11; i++)
    {
        p.posx.push_back(0);
        p.posy.push_back(0);
        p.posz.push_back(-(10-i));
    }

    // Magnetic rotation constants
    double qmdt2 = (p.q / p.m) * (dt / 2.0);
    double tx = qmdt2 * Bx;
    double t_mag2 = tx * tx;
    double sx = 2 * tx / (1 + t_mag2);

    for (int i = 0; i < steps; i++)
    {
        // Boris rotation in y-z plane
        double vpy = p.vy + p.vz * tx;
        double vpz = p.vz - p.vy * tx;

        p.vy += vpz * sx;
        p.vz -= vpy * sx;

        // Position update (x doesn't change)
        p.y += p.vy * dt;
        p.z += p.vz * dt;

        // Save trajectory
        if(std::fmod(i, save_step) == 0)
        {
            p.posx.push_back(p.x/mm); // stays constant
            p.posy.push_back(p.y/mm);
            p.posz.push_back(p.z/mm);
        }

        if (std::abs(p.y) >= magnet_height / 2 || (p.z/mm)<-0.05)
        {
            break;
        }
    }

    // Final speed and energy
    p.v = std::sqrt(p.vy * p.vy + p.vz * p.vz); // vx = 0
    p.energy = 0.5 * p.m * p.v * p.v;
}


void print_1d_vector(const std::vector<double>& v) {
    for (const auto& elem : v) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}


std::vector<double> concatenate(const std::initializer_list<std::vector<double>> &vectors)
{
    std::vector<double> result;

    // Reserve total size to avoid reallocations
    size_t total_size = 0;
    for (const auto &vec : vectors)
        total_size += vec.size();
    result.reserve(total_size);

    // Insert all elements
    for (const auto &vec : vectors)
        result.insert(result.end(), vec.begin(), vec.end());

    return result;
}


std::vector<double> linspace(double start, double end, int num) 
{
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result[i] = start + step * i;
    }
    return result;
}


std::pair<double, double> solve2x2(double a11, double a12, double a21, double a22, double b1, double b2) {
    double det = a11 * a22 - a12 * a21;
    if (std::abs(det) < 1e-12) throw std::runtime_error("Singular matrix");

    double alpha = (b1 * a22 - b2 * a12) / det;
    double beta  = (a11 * b2 - a21 * b1) / det;

    return {alpha, beta};
}


std::vector<double> generate_scaled_energy(double low_energy, double high_energy, int no_of_particles, double steep)
{
    std::vector<double> raw = linspace(0.0, 10.0, no_of_particles);
    for (double& val : raw) {
        val = std::exp(-val / steep);
    }

    double last_raw = raw.back();

    // Linear system to get alpha, beta
    auto [alpha, beta] = solve2x2(
        low_energy, (high_energy - low_energy),
        low_energy, last_raw * (high_energy - low_energy),
        high_energy, low_energy
    );

    std::vector<double> Energy(no_of_particles);
    for (int i = 0; i < no_of_particles; ++i) {
        Energy[i] = (alpha * low_energy + beta * raw[i] * (high_energy - low_energy)) * kev_to_joule;
    }

    return Energy;
}


void plotScatterPlot(const std::vector<double>& x, const std::vector<double>& y, const std::string& color, double pointSize) {
    if (x.size() != y.size()) {
        std::cerr << "Error: The vectors must have the same size!" << std::endl;
        return;
    }

    // Step 1: Write the data to a file
    std::ofstream dataFile("scatter_data.dat");
    for (size_t i = 0; i < x.size(); ++i) {
        dataFile << x[i] << " " << y[i] << "\n";
    }
    dataFile.close();

    // Step 2: Determine the plot limits with padding
    double xMin = *std::min_element(x.begin(), x.end());
    double xMax = *std::max_element(x.begin(), x.end());
    double yMin = *std::min_element(y.begin(), y.end());
    double yMax = *std::max_element(y.begin(), y.end());

    // Add padding to the axis limits (10% of the range)
    double xPadding = (xMax - xMin) * 0.1;
    double yPadding = (yMax - yMin) * 0.1;

    xMin -= xPadding;
    xMax += xPadding;
    yMin -= yPadding;
    yMax += yPadding;

    // Step 3: Open a pipe to GNUPLOT
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        // Send plotting commands to GNUPLOT
        fprintf(gnuplotPipe, "set title 'Scatter Plot'\n");
        fprintf(gnuplotPipe, "set xlabel 'X Axis'\n");
        fprintf(gnuplotPipe, "set ylabel 'Y Axis'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "set xrange [%f:%f]\n", xMin, xMax);  // Set x-axis limits with padding
        fprintf(gnuplotPipe, "set yrange [%f:%f]\n", yMin, yMax);  // Set y-axis limits with padding
        fprintf(gnuplotPipe, "plot 'scatter_data.dat' with points pt 7 ps %f lc rgb '%s'\n", pointSize, color.c_str());
        fflush(gnuplotPipe);

        // Wait for the user to close the plot window
        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();

        // Close the pipe to GNUPLOT
        pclose(gnuplotPipe);

        // Delete the data file
        if (remove("scatter_data.dat") != 0) {
            std::cerr << "Error deleting scatter_data.dat" << std::endl;
        } else {
            std::cout << "Deleted scatter_data.dat" << std::endl;
        }
    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }
}


void plotScatterPlot_multicolor(const std::vector<double>& x, const std::vector<double>& y, double pointSize, int groupSize) {
    if (x.size() != y.size()) {
        std::cerr << "Error: The vectors must have the same size!" << std::endl;
        return;
    }

    // Define a list of colors to cycle through
    std::vector<std::string> colors = {
        "red", "blue", "green", "magenta", "cyan", "orange", "black", "purple", "brown"
    };

    size_t totalPoints = x.size();
    size_t numGroups = (totalPoints + groupSize - 1) / groupSize;

    // Step 1: Write grouped data files
    for (size_t g = 0; g < numGroups; ++g) {
        std::ofstream dataFile("scatter_data_" + std::to_string(g) + ".dat");
        size_t start = g * groupSize;
        size_t end = std::min(start + groupSize, totalPoints);
        for (size_t i = start; i < end; ++i) {
            dataFile << x[i] << " " << y[i] << "\n";
        }
        dataFile.close();
    }

    // Step 2: Determine padded plot limits
    double xMin = *std::min_element(x.begin(), x.end());
    double xMax = *std::max_element(x.begin(), x.end());
    double yMin = *std::min_element(y.begin(), y.end());
    double yMax = *std::max_element(y.begin(), y.end());
    double xPadding = (xMax - xMin) * 0.1;
    double yPadding = (yMax - yMin) * 0.1;
    xMin -= xPadding; xMax += xPadding;
    yMin -= yPadding; yMax += yPadding;

    // Step 3: Plot with GNUPLOT
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title 'Scatter Plot with Color Groups'\n");
        fprintf(gnuplotPipe, "set xlabel 'Electric shift (mm)'\n");
        fprintf(gnuplotPipe, "set ylabel 'Magnetic shift (mm)'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "set xrange [%f:%f]\n", xMin, xMax);
        fprintf(gnuplotPipe, "set yrange [%f:%f]\n", yMin, yMax);
        fprintf(gnuplotPipe, "plot ");

        for (size_t g = 0; g < numGroups; ++g) {
            std::string color = colors[g % colors.size()];
            fprintf(gnuplotPipe, "'scatter_data_%zu.dat' with points pt 7 ps %f lc rgb '%s' title 'Group %zu'", 
                g, pointSize, color.c_str(), g + 1);
            if (g != numGroups - 1) {
                fprintf(gnuplotPipe, ", ");
            } else {
                fprintf(gnuplotPipe, "\n");
            }
        }

        fflush(gnuplotPipe);
        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();
        pclose(gnuplotPipe);

        // Delete all group files
        for (size_t g = 0; g < numGroups; ++g) {
            std::string filename = "scatter_data_" + std::to_string(g) + ".dat";
            if (remove(filename.c_str()) != 0) {
                std::cerr << "Error deleting " << filename << std::endl;
            }
        }

    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }
}


std::string jet_color(double t) 
{
    double r = std::clamp(1.5 - std::abs(4 * t - 3), 0.0, 1.0);
    double g = std::clamp(1.5 - std::abs(4 * t - 2), 0.0, 1.0);
    double b = std::clamp(1.5 - std::abs(4 * t - 1), 0.0, 1.0);

    auto toHex = [](double val) {
        int i = static_cast<int>(val * 255);
        char hex[3];
        snprintf(hex, sizeof(hex), "%02x", i);
        return std::string(hex);
    };

    return "#" + toHex(r) + toHex(g) + toHex(b);
}


void plotScatterPlot_multicolor_custom(const std::vector<double>& x, const std::vector<double>& y, double pointSize, const std::vector<size_t>& segmentSizes) 
{
    if (x.size() != y.size()) {
        std::cerr << "Error: The vectors must have the same size!" << std::endl;
        return;
    }

    // Step 1: Write segmented data files
    size_t currentIndex = 0;
    size_t groupCount = segmentSizes.size();
    for (size_t g = 0; g < groupCount; ++g) {
        size_t groupSize = segmentSizes[g];
        std::ofstream dataFile("scatter_data_" + std::to_string(g) + ".dat");
        for (size_t i = 0; i < groupSize && currentIndex < x.size(); ++i, ++currentIndex) {
            dataFile << x[currentIndex] << " " << y[currentIndex] << "\n";
        }
        dataFile.close();
    }

    // Step 2: Determine plot limits with padding
    double xMin = *std::min_element(x.begin(), x.end());
    double xMax = *std::max_element(x.begin(), x.end());
    double yMin = *std::min_element(y.begin(), y.end());
    double yMax = *std::max_element(y.begin(), y.end());
    double xPadding = (xMax - xMin) * 0.1;
    double yPadding = (yMax - yMin) * 0.1;
    xMin -= xPadding; xMax += xPadding;
    yMin -= yPadding; yMax += yPadding;

    // Step 3: Plot using Gnuplot with dynamic colormap
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title 'Scatter Plot with Jet Colormap'\n");
        fprintf(gnuplotPipe, "set xlabel 'Electric shift (mm)'\n");
        fprintf(gnuplotPipe, "set ylabel 'Magnetic shift (mm)'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "set xrange [%f:%f]\n", xMin, xMax);
        fprintf(gnuplotPipe, "set yrange [%f:%f]\n", yMin, yMax);
        fprintf(gnuplotPipe, "set size ratio -1\n");  // Ensure equal aspect ratio
        fprintf(gnuplotPipe, "plot ");

        for (size_t g = 0; g < groupCount; ++g) {
            double t = static_cast<double>(g) / std::max(groupCount - 1, size_t(1));
            std::string color = jet_color(t);
            fprintf(gnuplotPipe, "'scatter_data_%zu.dat' with points pt 7 ps %f lc rgb '%s' notitle",
                    g, pointSize, color.c_str());
            if (g != groupCount - 1) {
                fprintf(gnuplotPipe, ", ");
            } else {
                fprintf(gnuplotPipe, "\n");
            }
        }

        fflush(gnuplotPipe);
        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();
        pclose(gnuplotPipe);

        // Delete all data files
        for (size_t g = 0; g < groupCount; ++g) {
            std::string filename = "scatter_data_" + std::to_string(g) + ".dat";
            if (remove(filename.c_str()) != 0)
                std::cerr << "Error deleting " << filename << std::endl;
        }
    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }
}


void plotLinePlot_multicolor_custom(const std::vector<double>& x, const std::vector<double>& y, double lineWidth, const std::vector<size_t>& segmentSizes) 
{
    if (x.size() != y.size()) {
        std::cerr << "Error: The vectors must have the same size!" << std::endl;
        return;
    }

    // Step 1: Write segmented data files
    size_t currentIndex = 0;
    size_t groupCount = segmentSizes.size();
    for (size_t g = 0; g < groupCount; ++g) {
        size_t groupSize = segmentSizes[g];
        std::ofstream dataFile("line_data_" + std::to_string(g) + ".dat");
        for (size_t i = 0; i < groupSize && currentIndex < x.size(); ++i, ++currentIndex) {
            dataFile << x[currentIndex] << " " << y[currentIndex] << "\n";
        }
        dataFile.close();
    }

    // Step 2: Determine plot limits with padding
    double xMin = *std::min_element(x.begin(), x.end());
    double xMax = *std::max_element(x.begin(), x.end());
    double yMin = *std::min_element(y.begin(), y.end());
    double yMax = *std::max_element(y.begin(), y.end());
    double xPadding = (xMax - xMin) * 0.1;
    double yPadding = (yMax - yMin) * 0.1;
    xMin -= xPadding; xMax += xPadding;
    yMin -= yPadding; yMax += yPadding;

    // Step 3: Plot using Gnuplot with dynamic colormap
    FILE* gnuplotPipe = popen("gnuplot -persistent", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set title 'particle trajectory in ESM (Bx: %.2f T)'\n", Bx);
        fprintf(gnuplotPipe, "set xlabel 'z (mm)'\n");
        fprintf(gnuplotPipe, "set ylabel 'y (mm)'\n");
        fprintf(gnuplotPipe, "set grid\n");
        fprintf(gnuplotPipe, "set xrange [%f:%f]\n", xMin, xMax);
        fprintf(gnuplotPipe, "set yrange [%f:%f]\n", yMin, yMax);
        fprintf(gnuplotPipe, "set size ratio -1\n");  // Ensure equal aspect ratio
        fprintf(gnuplotPipe, "plot ");

        for (size_t g = 0; g < groupCount; ++g) {
            double t = static_cast<double>(g) / std::max(groupCount - 1, size_t(1));
            std::string color = jet_color(t);
            fprintf(gnuplotPipe, "'line_data_%zu.dat' with lines lw %f lc rgb '%s' notitle",
                    g, lineWidth, color.c_str());
            if (g != groupCount - 1) {
                fprintf(gnuplotPipe, ", ");
            } else {
                fprintf(gnuplotPipe, "\n");
            }
        }

        fflush(gnuplotPipe);
        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();
        pclose(gnuplotPipe);

        // Delete all data files
        for (size_t g = 0; g < groupCount; ++g) {
            std::string filename = "line_data_" + std::to_string(g) + ".dat";
            if (remove(filename.c_str()) != 0)
                std::cerr << "Error deleting " << filename << std::endl;
        }
    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }
}


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
    // gp_script << "set view equal xyz\n";  // Enforce equal aspect ratio
    gp_script << "set terminal wxt size 800,600\n";
    gp_script << "splot ";

    // Add each segment to the splot command with different colors
    for (size_t s = 0; s < num_segments; ++s) {
        std::string filename = "segment_" + std::to_string(s) + ".dat";
        // gp_script << "'" << filename << "' using 1:2:3 with points pointtype 7 pointsize 1.5 title 'Segment " << s << "'";
        gp_script << "'" << filename << "' using 1:2:3 with points pointtype 7 pointsize 1 notitle";

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
    double low_energy=100;
    double high_energy=20100;
    double steep = 2;
    double no_of_particle = 50;
    double t_res = 50000;

    std::vector<double> energy_array = generate_scaled_energy(low_energy, high_energy, no_of_particle, steep);
    std::vector<double> xs={},ys={},zs={};
    std::vector<size_t> segmentSizes = {};

    for(int i=0;i<no_of_particle;i++)
    {
        double x=0,y=0,z=0;
        double vx=0,vy=0;//,vz=sqrt(2*energy_array[i]/me);
        
        
        double gamma = 1+energy_array[i]/(me*c*c);
        double vz = c*sqrt(1-1/(gamma*gamma)); 
        double m = me*gamma;

        double t_max = 50*magnet_height/(vz);
        double dt = t_max/t_res;

        Particle p = Particle::__init__{.x=x, .y=y, .z=z, .vx=vx, .vy=vy, .vz=vz, .q=qe, .m=m}.__build__();
        propagator(p,t_max,dt,10);

        xs=concatenate({xs,p.posx});
        ys=concatenate({ys,p.posy});
        zs=concatenate({zs,p.posz});

        segmentSizes.push_back(p.posx.size());
        cout<<"progress: "<<(i+1)*100/no_of_particle<<"%"<<endl;
    }

    // plot_3d_scatter_segments(xs,ys,zs,static_cast<int>(xs.size()/no_of_particle));
    // plotScatterPlot_multicolor_custom(ys,zs,1,segmentSizes);
    plotLinePlot_multicolor_custom(zs,ys,1.5,segmentSizes);
    
    return 0;
}


