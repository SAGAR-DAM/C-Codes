#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

// Constants
constexpr double keV_to_Joule = 1;

// Linear interpolation like np.linspace

std::vector<double> linspace(double start, double end, int num) 
{
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result[i] = start + step * i;
    }
    return result;
}

// Solve 2x2 linear system: A * [alpha, beta] = b
std::pair<double, double> solve2x2(double a11, double a12, double a21, double a22, double b1, double b2) {
    double det = a11 * a22 - a12 * a21;
    if (std::abs(det) < 1e-12) throw std::runtime_error("Singular matrix");

    double alpha = (b1 * a22 - b2 * a12) / det;
    double beta  = (a11 * b2 - a21 * b1) / det;

    return {alpha, beta};
}

// Main energy generation function
std::vector<double> generate_scaled_energy(
    double low_energy,
    double high_energy,
    int no_of_particles,
    double steep)
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
        Energy[i] = (alpha * low_energy + beta * raw[i] * (high_energy - low_energy)) * keV_to_Joule;
    }

    return Energy;
}

// Main function — clean and minimal
int main() {
    double low_energy = 500.0;
    double high_energy = 40000.0;
    int no_of_particles = 1000;
    double steep = 2.0;

    std::vector<double> Energy = generate_scaled_energy(low_energy, high_energy, no_of_particles, steep);

    // Optional: print a few
    for (int i = 0; i < Energy.size(); ++i) {
        std::cout << "Energy[" << i << "] = " << Energy[i] << " J" << std::endl;
    }

    return 0;
}
