#include <iostream>
#include <vector>
#include <cmath>
#include <cstdint>  // Include this for int64_t and uint64_t
#include "D:/C++modules/matplotlibcpp/matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;

// Constants
const double m1 = 5.0;  // mass of the central body
const double m2 = 1.0;  // mass of the moving object
const double mu = (m1 * m2) / (m1 + m2);  // reduced mass
const double h = 1.0;   // angular momentum (constant)

// Define the central force F(r)
double F(double r) {
    double k = 1.0;  // strength of the force
    return -k / (1.0 + std::log(r));  // example log force
}

// Define the system of differential equations
void equations(double r, double r_dot, double theta, double &r_ddot, double &theta_dot) {
    r_ddot = F(r) / mu + (h * h) / (r * r * r);  // Radial acceleration
    theta_dot = h / (r * r);  // Angular velocity
}

// Runge-Kutta 4th order method (RK4)
void runge_kutta_4(double &r, double &r_dot, double &theta, double dt) {
    double k1_r, k1_r_dot, k1_theta;
    double k2_r, k2_r_dot, k2_theta;
    double k3_r, k3_r_dot, k3_theta;
    double k4_r, k4_r_dot, k4_theta;
    double r_ddot, theta_dot;

    // Step 1
    equations(r, r_dot, theta, r_ddot, theta_dot);
    k1_r = r_dot;
    k1_r_dot = r_ddot;
    k1_theta = theta_dot;

    // Step 2
    equations(r + 0.5 * k1_r * dt, r_dot + 0.5 * k1_r_dot * dt, theta + 0.5 * k1_theta * dt, r_ddot, theta_dot);
    k2_r = r_dot + 0.5 * k1_r_dot * dt;
    k2_r_dot = r_ddot;
    k2_theta = theta_dot;

    // Step 3
    equations(r + 0.5 * k2_r * dt, r_dot + 0.5 * k2_r_dot * dt, theta + 0.5 * k2_theta * dt, r_ddot, theta_dot);
    k3_r = r_dot + 0.5 * k2_r_dot * dt;
    k3_r_dot = r_ddot;
    k3_theta = theta_dot;

    // Step 4
    equations(r + k3_r * dt, r_dot + k3_r_dot * dt, theta + k3_theta * dt, r_ddot, theta_dot);
    k4_r = r_dot + k3_r_dot * dt;
    k4_r_dot = r_ddot;
    k4_theta = theta_dot;

    // Combine increments
    r += (dt / 6.0) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r);
    r_dot += (dt / 6.0) * (k1_r_dot + 2 * k2_r_dot + 2 * k3_r_dot + k4_r_dot);
    theta += (dt / 6.0) * (k1_theta + 2 * k2_theta + 2 * k3_theta + k4_theta);
}

int main() {
    // Initial conditions
    double r = 2.0;         // initial radial distance
    double r_dot = 0.2;     // initial radial velocity
    double theta = 0.0;     // initial angle
    double t_max = 400.0;   // total simulation time
    double dt = 0.1;        // time step
    int n_steps = t_max / dt;  // number of time steps

    // Vectors to store the trajectory
    std::vector<double> x_vals, y_vals;

    // Simulate using RK4
    for (int i = 0; i < n_steps; ++i) {
        // Convert polar coordinates (r, theta) to Cartesian (x, y)
        double x = r * cos(theta);
        double y = r * sin(theta);
        
        // Store the current position
        x_vals.push_back(x);
        y_vals.push_back(y);

        // Perform one RK4 step
        runge_kutta_4(r, r_dot, theta, dt);
    }

    // Set up the figure
    plt::figure();
    plt::xlim(-3, 3);
    plt::ylim(-3, 3);
    plt::plot(x_vals, y_vals, "b-");  // Trajectory of the moving object
    plt::plot({0}, {0}, "ro");  // Central body as a red dot

    // Show the plot
    plt::show();

    return 0;
}
