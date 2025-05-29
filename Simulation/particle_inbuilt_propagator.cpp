#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double Ex=12000;
double Ey=10000;
double Ez=0;
double Bx=1;
double By=1;
double Bz=0;

double cm=1e-2;
double mm=1e-3;
double qe=1.6e-19;
double mH=1.67e-27;

double kev_to_joule=1.6e-16;

struct Particle
{
    /* data */
    double x,y,z,vx,vy,vz,q,m,v,energy;
    Particle(   double x=0.0,
                double y=0.0,
                double z=0.0,
                double vx=0.0,
                double vy=0.0,
                double vz=0.0,
                double q=0.0,
                double m=1.0)
    {
        this->x=x;
        this->y=y;
        this->z=z;
        this->vx=vx;
        this->vy=vy;
        this->vz=vz;
        this->q=q;
        this->m=m;
        v = sqrt(vx*vx + vy*vy + vz*vz);
        energy = 0.5 * m * v * v;
    }

    struct __init__
    {
        double x=0.0;
        double y=0.0;
        double z=0.0;
        double vx=0.0;
        double vy=0.0;
        double vz=0.0;
        double q=0.0;
        double m=1.0;

        Particle __build__() const 
        { 
            return Particle(x,y,z,vx,vy,vz,q,m); 
        }
    };


    std::vector<double> posx;
    std::vector<double> posy;
    std::vector<double> posz;


    void propagator(double t_max = 0, double dt = 0.001)
    {
        double t = 0.0;
        int steps = static_cast<int>(t_max / dt);

        // Constants for field
        double qmdt2 = (q / m) * (dt / 2.0);

        for (int i = 0; i < steps; i++)
        {
            // Half acceleration from E
            double vx_minus = vx + qmdt2 * Ex;
            double vy_minus = vy + qmdt2 * Ey;
            double vz_minus = vz + qmdt2 * Ez;

            // t vector
            double tx = qmdt2 * Bx;
            double ty = qmdt2 * By;
            double tz = qmdt2 * Bz;

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
            vx = vx_plus + qmdt2 * Ex;
            vy = vy_plus + qmdt2 * Ey;
            vz = vz_plus + qmdt2 * Ez;

            // Position update using full-step velocity
            x += vx * dt;
            y += vy * dt;
            z += vz * dt;

            // Save position
            posx.push_back(x);
            posy.push_back(y);
            posz.push_back(z);
        }

        // Final velocity and energy
        v = sqrt(vx * vx + vy * vy + vz * vz);
        energy = 0.5 * m * v * v;
    }

};


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
    double energy = 1 * kev_to_joule;
    double vx = 0, vy = 0, vz = pow(2 * energy / mH, 0.5);

    cout << "before assign: " << vx << "  " << vy << "  " << vz << endl;
    Particle p = Particle::__init__{.vx = vx, .vy = vy, .vz = vz, .q = qe, .m = mH}.__build__();
    cout << "after assign: " << p.vx << "  " << p.vy << "  " << p.vz << endl;
    cout << "p.v: " << p.v << endl;
    double t_max = 100 * cm / p.v;
    int res_t = 20000;
    double dt = t_max / res_t;

    cout << "t: " << t_max << " dt: " << dt << endl;

    p.propagator(t_max, dt);
    cout << p.vx << "   " << p.vy << "   " << p.vz << endl;
    cout<<"p.energy: "<<p.energy/kev_to_joule<<endl;
    // print_1d_vector(p.posx);

    // Generate time vector
    // std::vector<double> t_vec(p.posx.size());
    // double t=0;
    // for (size_t i = 0; i < t_vec.size(); ++i, t += dt) {
    //     t_vec[i] = t;
    // }

    // cout<<"p.m: "<<p.m<<" p.q: "<<p.q<<" p.energy: "<<p.energy/kev_to_joule<<" p.v: "<<p.v<<endl;

    std::vector<double> t_vec(p.posx.size());
    double t=0;
    for (size_t i = 0; i < t_vec.size(); ++i, t += dt) {
        t_vec[i] = t;
    }
    
    plot_3d_line(p.posx,p.posy,p.posz);

    
    return 0;
}
