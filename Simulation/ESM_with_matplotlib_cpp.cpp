// # -*- coding: utf-8 -*-
// """
// Created on Thu May 29 19:40:44 2025

// @author: mrsag
// """

#include "C:\\Users\\mrsag\\AppData\\Local\\Programs\\Python\\Python311\\include\\Python.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <limits>

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

double magnet_height = 20 * mm;

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

void propagator(Particle &p, double t_max = 0, double dt = 0.001, int save_step = 50)
{
    int steps = static_cast<int>(t_max / dt);
    for (int i = 0; i < 11; i++)
    {
        p.posx.push_back(0);
        p.posy.push_back(0);
        p.posz.push_back(-(10 - i));
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
        if (std::fmod(i, save_step) == 0)
        {
            p.posx.push_back(p.x / mm); // stays constant
            p.posy.push_back(p.y / mm);
            p.posz.push_back(p.z / mm);
        }

        if (std::abs(p.y) >= magnet_height / 2 || (p.z / mm) < -0.05)
        {
            break;
        }
    }

    // Final speed and energy
    p.v = std::sqrt(p.vy * p.vy + p.vz * p.vz); // vx = 0
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
    for (int i = 0; i < num; ++i)
    {
        result[i] = start + step * i;
    }
    return result;
}

std::pair<double, double> solve2x2(double a11, double a12, double a21, double a22, double b1, double b2)
{
    double det = a11 * a22 - a12 * a21;
    if (std::abs(det) < 1e-12)
        throw std::runtime_error("Singular matrix");

    double alpha = (b1 * a22 - b2 * a12) / det;
    double beta = (a11 * b2 - a21 * b1) / det;

    return {alpha, beta};
}

std::vector<double> generate_scaled_energy(double low_energy, double high_energy, int no_of_particles, double steep)
{
    std::vector<double> raw = linspace(0.0, 10.0, no_of_particles);
    for (double &val : raw)
    {
        val = std::exp(-val / steep);
    }

    double last_raw = raw.back();

    // Linear system to get alpha, beta
    auto [alpha, beta] = solve2x2(
        low_energy, (high_energy - low_energy),
        low_energy, last_raw * (high_energy - low_energy),
        high_energy, low_energy);

    std::vector<double> Energy(no_of_particles);
    for (int i = 0; i < no_of_particles; ++i)
    {
        Energy[i] = (alpha * low_energy + beta * raw[i] * (high_energy - low_energy)) * kev_to_joule;
    }

    return Energy;
}

// Converts one vector and injects it into the __main__ module under a given name
void convert_and_inject_vector(const std::vector<double> &vec, const char *name)
{
    PyObject *pyList = PyList_New(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
    {
        PyList_SetItem(pyList, i, PyFloat_FromDouble(vec[i]));
    }

    PyObject *mainModule = PyImport_AddModule("__main__");
    PyObject *mainDict = PyModule_GetDict(mainModule);
    PyDict_SetItemString(mainDict, name, pyList);

    Py_DECREF(pyList);
}

// Helper function to expand vector-name pairs by index
template <typename TupleVecs, typename TupleNames, std::size_t... Is>
void convert_and_inject_vectors_impl(const TupleVecs &vecs, const TupleNames &names, std::index_sequence<Is...>)
{
    (convert_and_inject_vector(std::get<Is>(vecs), std::get<Is>(names)), ...);
}

// Main function: accepts tuples of vectors and names
template <typename... Vectors, typename... Names>
void convert_and_inject_vectors(const std::tuple<Vectors...> &vecs, const std::tuple<Names...> &names)
{
    static_assert(sizeof...(Vectors) == sizeof...(Names), "Number of vectors and names must match.");
    convert_and_inject_vectors_impl(vecs, names, std::index_sequence_for<Vectors...>{});
}





int main()
{

    Py_Initialize();
    PyRun_SimpleString("import matplotlib.pyplot as plt\nimport numpy as np");


    double low_energy=100;
    double high_energy=20100;
    double steep = 2;
    double no_of_particle = 50;
    double t_res = 50000;

    std::vector<double> energy_array = generate_scaled_energy(low_energy, high_energy, no_of_particle, steep);

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

        convert_and_inject_vectors(std::make_tuple(p.posz,p.posy),
                               std::make_tuple("z","y"));

        const char* pyCode = R"(
plt.plot(z, y)
        )";
        PyRun_SimpleString(pyCode);
        cout<<"progress: "<<(i+1)*100/no_of_particle<<"%"<<endl;
    }

    const char* pyCode1 = R"(
plt.grid(True)
plt.gca().set_aspect('equal')
plt.show()
        )";    

    PyRun_SimpleString(pyCode1);

    Py_Finalize();
    return 0;

}
