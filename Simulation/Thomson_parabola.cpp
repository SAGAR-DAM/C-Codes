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

double Ex = 10000/(1.3*cm);
double Ey = 0;
double Ez = 0;
double Bx = 0.45;
double By = 0;
double Bz = 0;

double V = 12000 ; // Potential difference in between the cathode and anode (V)   (in Volts)
double d = 1.3*cm ; // distance over which the electric field is being applied (in cm)
double E_field = V/d ; // electric field (in Volts/m))
double B_field = 0.45 ; // Magnetic Field (B)  (in Tesla)
double lE = 12*cm ; // Electric field length (lE)  (in cm)
double lB = 5*cm  ; // Magnetic field length (lB)  (in cm)
double DE = 18.8*cm ; // Distance of screen from electric field region endpoint (DE)  (in cm)
double DB = (32.3+1.75)*cm ; // Distance of screen from magnetic field region endpoint (DB)  (in cm)


struct Particle
{
    /* data */
    double x, y, z, vx, vy, vz, q, m;
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



    double v = sqrt(vx*vx+vy*vy+vz*vz);
    double energy = 0.5*m*v*v;


    void propagator(double t_max = 0, double dt = 0.001)
    {
        double t = 0.0;
        int steps = static_cast<int>(t_max / dt);
        double Exf=0,Eyf=0,Ezf=0,Bxf=0,Byf=0,Bzf=0;
        bool broken=false;

        // Constants for field
        double qmdt2 = (q / m) * (dt / 2.0);

        for (int i = 0; i < steps; i++)
        {   
            if(abs(z)>=DE && abs(z)<=lE+DE)
            {
                Exf=Ex;
                Eyf=Ey;
                Ezf=Ez;
            }
            else
            {
                Exf=0;
                Eyf=0;
                Ezf=0;
            }
            if(abs(z)>=DB && abs(z)<=lB+DB)
            {
                Bxf=Bx;
                Byf=By;
                Bzf=Bz;
            }
            else
            {
                Bxf=0;
                Byf=0;
                Bzf=0;
            }

            // Half acceleration from E
            double vx_minus = vx + qmdt2 * Exf;
            double vy_minus = vy + qmdt2 * Eyf;
            double vz_minus = vz + qmdt2 * Ezf;

            // t vector
            double tx = qmdt2 * Bxf;
            double ty = qmdt2 * Byf;
            double tz = qmdt2 * Bzf;

            // t^2
            double t_mag2 = tx * tx + ty * ty + tz * tz;

            // s vector
            double sx = 2 * tx / (1 + t_mag2);
            double sy = 2 * ty / (1 + t_mag2);
            double sz = 2 * tz / (1 + t_mag2);

            // v' = v_minus + v_minus x t
            double vpx = vx_minus + (vy_minus * tz - vz_minus * ty);
            double vpy = vy_minus + (vz_minus * tx - vx_minus * tz);
            double vpz = vz_minus + (vx_minus * ty - vy_minus * tx);

            // v_plus = v_minus + v' x s
            double vx_plus = vx_minus + (vpy * sz - vpz * sy);
            double vy_plus = vy_minus + (vpz * sx - vpx * sz);
            double vz_plus = vz_minus + (vpx * sy - vpy * sx);

            // Final velocity (after second half E kick)
            vx = vx_plus + qmdt2 * Exf;
            vy = vy_plus + qmdt2 * Eyf;
            vz = vz_plus + qmdt2 * Ezf;

            // Position update using full-step velocity
            x += vx * dt;
            y += vy * dt;
            z += vz * dt;

            v = sqrt(vx * vx + vy * vy + vz * vz);
            energy = 0.5 * m * v * v;

            if(z>=0)
            {
                broken=true;
                break;
            }

        }

    }
    

};



struct Proton : public Particle
{
    Proton(double x = 0.0,
           double y = 0.0,
           double z = 0.0,
           double vx = 0.0,
           double vy = 0.0,
           double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, qe, mH)
    {
    }
};


struct C1 : public Particle
{
    C1(double x = 0.0,
           double y = 0.0,
           double z = 0.0,
           double vx = 0.0,
           double vy = 0.0,
           double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, qe, 12*mH)
    {
    }
};


struct C2 : public Particle
{
    C2(double x = 0.0,
           double y = 0.0,
           double z = 0.0,
           double vx = 0.0,
           double vy = 0.0,
           double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 2*qe, 12*mH)
    {
    }
};


struct C3 : public Particle
{
    C3(double x = 0.0,
           double y = 0.0,
           double z = 0.0,
           double vx = 0.0,
           double vy = 0.0,
           double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 3*qe, 12*mH)
    {
    }
};


struct C4 : public Particle
{
    C4(double x = 0.0,
           double y = 0.0,
           double z = 0.0,
           double vx = 0.0,
           double vy = 0.0,
           double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 4*qe, 12*mH)
    {
    }
};


struct C5 : public Particle
{
    C5(double x = 0.0,
           double y = 0.0,
           double z = 0.0,
           double vx = 0.0,
           double vy = 0.0,
           double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 5*qe, 12*mH)
    {
    }
};


struct C6 : public Particle
{
    C6(double x = 0.0,
           double y = 0.0,
           double z = 0.0,
           double vx = 0.0,
           double vy = 0.0,
           double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 6*qe, 12*mH)
    {
    }
};


struct O1 : public Particle {
    O1(double x = 0.0, double y = 0.0, double z = 0.0,
       double vx = 0.0, double vy = 0.0, double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 1 * qe, 16 * mH) {}
};


struct O2 : public Particle {
    O2(double x = 0.0, double y = 0.0, double z = 0.0,
       double vx = 0.0, double vy = 0.0, double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 2 * qe, 16 * mH) {}
};


struct O3 : public Particle {
    O3(double x = 0.0, double y = 0.0, double z = 0.0,
       double vx = 0.0, double vy = 0.0, double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 3 * qe, 16 * mH) {}
};


struct O4 : public Particle {
    O4(double x = 0.0, double y = 0.0, double z = 0.0,
       double vx = 0.0, double vy = 0.0, double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 4 * qe, 16 * mH) {}
};


struct O5 : public Particle {
    O5(double x = 0.0, double y = 0.0, double z = 0.0,
       double vx = 0.0, double vy = 0.0, double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 5 * qe, 16 * mH) {}
};


struct O6 : public Particle {
    O6(double x = 0.0, double y = 0.0, double z = 0.0,
       double vx = 0.0, double vy = 0.0, double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 6 * qe, 16 * mH) {}
};


struct O7 : public Particle {
    O7(double x = 0.0, double y = 0.0, double z = 0.0,
       double vx = 0.0, double vy = 0.0, double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 7 * qe, 16 * mH) {}
};


struct O8 : public Particle {
    O8(double x = 0.0, double y = 0.0, double z = 0.0,
       double vx = 0.0, double vy = 0.0, double vz = 0.0)
        : Particle(x, y, z, vx, vy, vz, 8 * qe, 16 * mH) {}
};


void print_1d_vector(const std::vector<double>& v) {
    for (const auto& elem : v) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
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



int main()
{
    std::vector<std::string> particle_names = {"proton", "C1", "C2", "C3", "C4", "C5", "C6","O1", "O2", "O3", "O4", "O5", "O6", "O7", "O8"};
    std::vector<double> low_energy_array  = {1200,8000,8000,8000,8000,8000,8000,10000,10000,10000,10000,10000,10000,10000,10000};
    std::vector<double> high_energy_array = {40000,100000,100000,100000,100000,100000,100000,180000,180000,180000,180000,180000,180000,180000,180000};
    std::vector<double> masses = {mH,12*mH,12*mH,12*mH,12*mH,12*mH,12*mH,16*mH,16*mH,16*mH,16*mH,16*mH,16*mH,16*mH,16*mH};
    std::vector<double> charges = {qe,1*qe,2*qe,3*qe,4*qe,5*qe,6*qe,1*qe,2*qe,3*qe,4*qe,5*qe,6*qe,7*qe,8*qe};
    double steep = 2.0;
    double no_of_particle = 5000;
    double t_res = 10000;

    std::vector<double> posx, posy; 

    double progress = 0;

    for(int i=0;i<particle_names.size();i++)
    {
        std::vector<double> energy_array = generate_scaled_energy(low_energy_array[i], high_energy_array[i], no_of_particle, steep);
        for(int j=0; j<energy_array.size();j++)
        {
            double x=0,y=0,z=std::max(lE + DE, lB + DB)*(-1);
            double vx=0,vy=0,vz=sqrt(2*energy_array[j]/masses[i]);
            double t_max = 5*std::max(lE + DE, lB + DB)/vz;
            double dt = t_max/t_res;

            Particle p = Particle::__init__{.x=x, .y=y, .z=z, .vx=vx, .vy=vy, .vz=vz, .q=charges[i], .m=masses[i]}.__build__();
            p.propagator(t_max,dt);
            posx.push_back(p.x/mm);
            posy.push_back(p.y/mm);

            if(std::fmod(progress / (no_of_particle * particle_names.size()) * 100, 1) == 0)
            {
                cout<<"progress: "<<progress/(no_of_particle*particle_names.size())*100<<" %"<<endl;
            }
            progress++;
        }

    }
    cout<<"progress: "<<100<<" %"<<endl;
    posx.push_back(0);
    posy.push_back(0);
    
    plotScatterPlot_multicolor(posx,posy,1,static_cast<int>(no_of_particle));



    return 0;
}