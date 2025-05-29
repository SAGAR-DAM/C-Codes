#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

double Ex = 12000;
double Ey = 10000;
double Ez = 0;
double Bx = 1;
double By = 1;
double Bz = 0;

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

// void propagator(Particle& p, double t_max=0, double dt=0.001)
// {
//     double Fx=0,Fy=0,Fz=0,ax=0,ay=0,az=0;
//     double t=0;
//     bool broke=false;

//     double t_steps = t_max/dt;
//     int timesteps = t_steps;

//     cout<<t_steps<<endl;

//     for(int i=0; i<timesteps; i++)
//     {
//         t += dt;

//         Fx = p.q*(Ex+(p.vy*Bz-p.vz*By));
//         Fy = p.q*(Ey+(p.vz*Bx-p.vx*Bz));
//         Fz = p.q*(Ez+(p.vx*By-p.vy*Bx));

//         ax = Fx / p.m;
//         ay = Fy / p.m;
//         az = Fz / p.m;

//         p.x += ((p.vx * dt) + (ax * dt*dt / 2));
//         p.y += ((p.vy * dt) + (ay * dt*dt / 2));
//         p.z += ((p.vz * dt) + (az * dt*dt / 2));

//         p.vx += (ax * dt);
//         p.vy += (ay * dt);
//         p.vz += (az * dt);

//         p.posx.push_back(p.x);
//         p.posy.push_back(p.y);
//         p.posz.push_back(p.z);

//         // cout<<"step: "<<i<<" ,time: "<<t/t_max<<endl;
//     }

//     p.v = sqrt(p.vx*p.vx + p.vy*p.vy + p.vz*p.vz);
//     p.energy = 0.5 * p.m * p.v * p.v;
// }

void propagator(Particle &p, double t_max = 0, double dt = 0.001)
{
    double t = 0.0;
    int steps = static_cast<int>(t_max / dt);

    // Constants for field
    double qmdt2 = (p.q / p.m) * (dt / 2.0);

    for (int i = 0; i < steps; i++)
    {
        // Half acceleration from E
        double vx_minus = p.vx + qmdt2 * Ex;
        double vy_minus = p.vy + qmdt2 * Ey;
        double vz_minus = p.vz + qmdt2 * Ez;

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
        p.vx = vx_plus + qmdt2 * Ex;
        p.vy = vy_plus + qmdt2 * Ey;
        p.vz = vz_plus + qmdt2 * Ez;

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

void print_1d_vector(const std::vector<double> &v)
{
    for (const auto &elem : v)
    {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
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

    propagator(p, t_max, dt);
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
    return 0;
}
