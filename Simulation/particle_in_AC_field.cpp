#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
// #include "D:\C++ codes\matplotlibcpp\include\matplotlibcpp.h"

// namespace plt = matplotlibcpp;


using namespace std;

double Ex = 0;
double Ey = 10000;
double Ez = 0;
double Bx = 0;
double By = 2;
double Bz = 0;

double omegax = 0e6;
double omegay = 0e6;
double omegaz = 0e6;

double cm = 1e-2;
double mm = 1e-3;
double qe = 1.60217663e-19;
double me = 9.10938e-31;
double mH = 1.67262192e-27;

double kev_to_joule = 1.660217663e-16;

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


void propagator(Particle &p, double t_max = 0, double dt = 0.001)
{
    double t = 0.0;
    int steps = static_cast<int>(t_max / dt);

    // Constants for field
    double qmdt2 = (p.q / p.m) * (dt / 2.0);
    double Exf=Ex,Eyf=Ey,Ezf=Ez,Bxf=Bx,Byf=By,Bzf=Bz;

    for (int i = 0; i < steps; i++)
    {
        Exf = Ex*cos(3.14159*i/steps);
        Eyf = Ey*cos(6*3.14159*i/steps);
        Ezf = Ex*cos(3.14159*i/steps);

        Byf = By*cos(2*3.14159*i/steps);

        // Half acceleration from E
        double vx_minus = p.vx + qmdt2 * Exf;
        double vy_minus = p.vy + qmdt2 * Eyf;
        double vz_minus = p.vz + qmdt2 * Ezf;

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
        p.vx = vx_plus + qmdt2 * Exf;
        p.vy = vy_plus + qmdt2 * Eyf;
        p.vz = vz_plus + qmdt2 * Ezf;

        // Position update using full-step velocity
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;

        // Save position
        p.posx.push_back(p.x);
        p.posy.push_back(p.y);
        p.posz.push_back(p.z);
    }

    // Final velocity and energy
    p.v = sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
    p.energy = 0.5 * p.m * p.v * p.v;
}


void print_1d_vector(const std::vector<double>& v) {
    for (const auto& elem : v) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}


void plot_3d_line(std::vector<double> posx, std::vector<double> posy, std::vector<double> posz)
{
    // Step 1: Write data to a file for Gnuplot
    std::ofstream data("3d_particle_data.dat");
    for (size_t i = 0; i < posx.size(); ++i) {
        // Writing x, y, z data points to the file
        data << posx[i] << " " << posy[i] << " " << posz[i] << "\n";
    }
    data.close();

    // Step 2: Create the Gnuplot script to plot 3D data
    std::ofstream gp_script("plot_3d.gp");
    gp_script << "set title '3D Particle Trajectory'\n";
    gp_script << "set xlabel 'X Position (m)'\n";
    gp_script << "set ylabel 'Y Position (m)'\n";
    gp_script << "set zlabel 'Z Position (m)'\n";
    gp_script << "set grid\n";
    gp_script << "set terminal wxt size 800,600\n";  // Interactive window
    gp_script << "splot '3d_particle_data.dat' using 1:2:3 with lines title 'Particle Path'\n";  // 3D line plot
    gp_script.close();

    // Step 3: Execute the Gnuplot script to generate the plot
    FILE *gnuplotPipe = popen("gnuplot", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "load 'plot_3d.gp'\n");
        fflush(gnuplotPipe);

        std::cout << "Press enter to exit..." << std::endl;
        std::cin.get();  // Wait for user input to keep the plot open

        pclose(gnuplotPipe);  // Close the pipe
    } else {
        std::cerr << "Error: Could not open pipe to GNUPLOT." << std::endl;
    }

    // Optionally, delete the data and script files after the plot is shown
    if (remove("3d_particle_data.dat") != 0) {
        std::cerr << "Error deleting 3d_particle_data.dat" << std::endl;
    } else {
        std::cout << "Deleted 3d_particle_data.dat" << std::endl;
    }

    if (remove("plot_3d.gp") != 0) {
        std::cerr << "Error deleting plot_3d.gp" << std::endl;
    } else {
        std::cout << "Deleted plot_3d.gp" << std::endl;
    }
}



int main()
{
    double energy = 1*kev_to_joule;
    double vx=0, vy=0, vz=pow(2*energy/mH,0.5);

    Particle p = Particle::__init__{.vx=vx,.vy=vy,.vz=vz,.q=qe,.m=mH}.__build__();
    double t_max = 100*cm/p.v;
    int res_t=10000;
    double dt=t_max/res_t;

    propagator(p,t_max,dt);

    // print_1d_vector(p.posx);

    // Generate time vector
    std::vector<double> t_vec(p.posx.size());
    double t=0;
    for (size_t i = 0; i < t_vec.size(); ++i, t += dt) {
        t_vec[i] = t;
    }
    
    plot_3d_line(p.posx,p.posy,p.posz);

    return 0;
}
