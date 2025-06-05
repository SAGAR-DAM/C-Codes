#include <iomanip>
#include <limits>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <iostream>
#include <algorithm>
#include "C:\\Users\\mrsag\\AppData\\Local\\Programs\\Python\\Python311\\include\\Python.h"

using namespace std;


class SimulationBox3D {
public:
    int nx, ny, nz;
    double lx, ly, lz;
    double dx, dy, dz;
    
    SimulationBox3D(int nx, int ny, int nz,
                    double lx, double ly, double lz,
                    double potential_offset = 0.0);

    void addSphere(double cx, double cy, double cz, double radius, double potential_value);
    void addBox(double x0, double y0, double z0,
            double x1, double y1, double z1, double potential_value);

    void addCylinder(double cx, double cy, double cz, double radius, double height, char axis, double potential_value);
    void addHollowPipe(double cx, double cy, double cz, double radius, double thickness, double height, char axis, double potential_value);
    void addEllipsoid(double cx, double cy, double cz, double rx, double ry, double rz, double potential_value);
    void addHyperboloid(double cx, double cy, double cz, double a, double b, double c, double waist, char axis, double potential_value);
    void addPlane(double A, double B, double C, double D, double thickness, double potential_value);



    void solve(int max_iter = 1000, double tol = 1e-4, const std::string& method = "jacobi");

    // Add geometry will be added next (e.g., addBox, addSphere)

    const std::vector<double>& getPotential() const { return potential; }


    std::vector<double> potential;
    std::vector<bool> fixed_mask;

    int index(int i, int j, int k) const {
        return i * ny * nz + j * nz + k;
    }



    void applyJacobi();
    void applyGaussSeidel();
};





SimulationBox3D::SimulationBox3D(int nx, int ny, int nz,
                                 double lx, double ly, double lz,
                                 double potential_offset)
    : nx(nx), ny(ny), nz(nz), lx(lx), ly(ly), lz(lz)
{
    dx = lx / (nx - 1);
    dy = ly / (ny - 1);
    dz = lz / (nz - 1);

    int total_size = nx * ny * nz;
    potential.resize(total_size, potential_offset);
    fixed_mask.resize(total_size, false);
}

void SimulationBox3D::solve(int max_iter, double tol, const std::string& method) {
    std::vector<double> V_old = potential;
    std::vector<double> V_new = potential;

    for (int iter = 0; iter < max_iter; ++iter) {
        V_old = potential;

        if (method == "jacobi") {
            applyJacobi();
        } else if (method == "gauss-seidel") {
            applyGaussSeidel();
        } else {
            throw std::runtime_error("Unknown method");
        }

        // Check for convergence
        double max_diff = 0.0;
        for (size_t i = 0; i < potential.size(); ++i) {
            max_diff = std::max(max_diff, std::abs(potential[i] - V_old[i]));
        }

        std::cout << "Iteration " << iter << ", Max Diff = " << max_diff << std::endl;

        if (max_diff < tol) {
            std::cout << "Converged at iteration " << iter << std::endl;
            break;
        }
    }
}



void SimulationBox3D::addSphere(double cx, double cy, double cz, double radius, double potential_value) {
    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        for (int j = 0; j < ny; ++j) {
            double y = j * dy;
            for (int k = 0; k < nz; ++k) {
                double z = k * dz;
                double dist = std::sqrt((x - cx)*(x - cx) +
                                        (y - cy)*(y - cy) +
                                        (z - cz)*(z - cz));
                if (dist <= radius) {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}

void SimulationBox3D::addBox(double x0, double y0, double z0,
                             double x1, double y1, double z1, double potential_value) {
    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        if (x < x0 || x > x1) continue;
        for (int j = 0; j < ny; ++j) {
            double y = j * dy;
            if (y < y0 || y > y1) continue;
            for (int k = 0; k < nz; ++k) {
                double z = k * dz;
                if (z < z0 || z > z1) continue;

                int idx = index(i, j, k);
                potential[idx] = potential_value;
                fixed_mask[idx] = true;
            }
        }
    }
}




void SimulationBox3D::addCylinder(double cx, double cy, double cz, double radius, double height, char axis, double potential_value) {
    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        for (int j = 0; j < ny; ++j) {
            double y = j * dy;
            for (int k = 0; k < nz; ++k) {
                double z = k * dz;

                bool inside = false;
                if (axis == 'z') {
                    double r2 = (x - cx)*(x - cx) + (y - cy)*(y - cy);
                    inside = (r2 <= radius*radius && z >= cz && z <= cz + height);
                } else if (axis == 'x') {
                    double r2 = (y - cy)*(y - cy) + (z - cz)*(z - cz);
                    inside = (r2 <= radius*radius && x >= cx && x <= cx + height);
                } else if (axis == 'y') {
                    double r2 = (x - cx)*(x - cx) + (z - cz)*(z - cz);
                    inside = (r2 <= radius*radius && y >= cy && y <= cy + height);
                }

                if (inside) {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}





void SimulationBox3D::addHollowPipe(double cx, double cy, double cz, double radius, double thickness, double height, char axis, double potential_value) {
    double r_outer2 = (radius + thickness / 2.0) * (radius + thickness / 2.0);
    double r_inner2 = (radius - thickness / 2.0) * (radius - thickness / 2.0);

    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        for (int j = 0; j < ny; ++j) {
            double y = j * dy;
            for (int k = 0; k < nz; ++k) {
                double z = k * dz;

                bool inside = false;
                if (axis == 'z') {
                    double r2 = (x - cx)*(x - cx) + (y - cy)*(y - cy);
                    inside = (r2 >= r_inner2 && r2 <= r_outer2 && z >= cz && z <= cz + height);
                } else if (axis == 'x') {
                    double r2 = (y - cy)*(y - cy) + (z - cz)*(z - cz);
                    inside = (r2 >= r_inner2 && r2 <= r_outer2 && x >= cx && x <= cx + height);
                } else if (axis == 'y') {
                    double r2 = (x - cx)*(x - cx) + (z - cz)*(z - cz);
                    inside = (r2 >= r_inner2 && r2 <= r_outer2 && y >= cy && y <= cy + height);
                }

                if (inside) {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}





void SimulationBox3D::addEllipsoid(double cx, double cy, double cz, double rx, double ry, double rz, double potential_value) {
    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        for (int j = 0; j < ny; ++j) {
            double y = j * dy;
            for (int k = 0; k < nz; ++k) {
                double z = k * dz;
                double value = ((x - cx) / rx) * ((x - cx) / rx)
                             + ((y - cy) / ry) * ((y - cy) / ry)
                             + ((z - cz) / rz) * ((z - cz) / rz);
                if (value <= 1.0) {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}





void SimulationBox3D::addHyperboloid(double cx, double cy, double cz, double a, double b, double c, double waist, char axis, double potential_value) {
    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        for (int j = 0; j < ny; ++j) {
            double y = j * dy;
            for (int k = 0; k < nz; ++k) {
                double z = k * dz;
                double val = 0.0;
                if (axis == 'x') {
                    val = -((x - cx)/a)*((x - cx)/a) + ((y - cy)/b)*((y - cy)/b) + ((z - cz)/c)*((z - cz)/c);
                } else if (axis == 'y') {
                    val = ((x - cx)/a)*((x - cx)/a) - ((y - cy)/b)*((y - cy)/b) + ((z - cz)/c)*((z - cz)/c);
                } else if (axis == 'z') {
                    val = ((x - cx)/a)*((x - cx)/a) + ((y - cy)/b)*((y - cy)/b) - ((z - cz)/c)*((z - cz)/c);
                }

                if (val <= waist * waist) {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}





void SimulationBox3D::addPlane(double A, double B, double C, double D, double thickness, double potential_value) {
    double norm = std::sqrt(A*A + B*B + C*C);
    for (int i = 0; i < nx; ++i) {
        double x = i * dx;
        for (int j = 0; j < ny; ++j) {
            double y = j * dy;
            for (int k = 0; k < nz; ++k) {
                double z = k * dz;
                double dist = (A * x + B * y + C * z + D) / norm;
                if (std::abs(dist) <= thickness / 2.0) {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}



void SimulationBox3D::applyJacobi()
{
    std::vector<double> new_potential = potential;

    // Interior points
    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            for (int k = 1; k < nz - 1; ++k)
            {
                int idx = index(i, j, k);
                if (!fixed_mask[idx])
                {
                    new_potential[idx] = (1.0 / 6.0) * (potential[index(i + 1, j, k)] + potential[index(i - 1, j, k)] +
                                                       potential[index(i, j + 1, k)] + potential[index(i, j - 1, k)] +
                                                       potential[index(i, j, k + 1)] + potential[index(i, j, k - 1)]);
                }
            }
        }
    }

    // Face points
    for (int j = 1; j < ny - 1; ++j)
    {
        for (int k = 1; k < nz - 1; ++k)
        {
            int idx = index(0, j, k);
            if (!fixed_mask[idx])
                new_potential[idx] = (1.0 / 5.0) * (potential[index(1, j, k)] + potential[index(0, j + 1, k)] +
                                                   potential[index(0, j - 1, k)] + potential[index(0, j, k + 1)] +
                                                   potential[index(0, j, k - 1)]);

            int idx1 = index(nx - 1, j, k);
            if (!fixed_mask[idx1])
                new_potential[idx1] = (1.0 / 5.0) * (potential[index(nx - 2, j, k)] + potential[index(nx - 1, j + 1, k)] +
                                                    potential[index(nx - 1, j - 1, k)] + potential[index(nx - 1, j, k + 1)] +
                                                    potential[index(nx - 1, j, k - 1)]);
        }
    }

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int k = 1; k < nz - 1; ++k)
        {
            int idx = index(i, 0, k);
            if (!fixed_mask[idx])
                new_potential[idx] = (1.0 / 5.0) * (potential[index(i + 1, 0, k)] + potential[index(i - 1, 0, k)] +
                                                   potential[index(i, 1, k)] + potential[index(i, 0, k + 1)] +
                                                   potential[index(i, 0, k - 1)]);

            int idx1 = index(i, ny - 1, k);
            if (!fixed_mask[idx1])
                new_potential[idx1] = (1.0 / 5.0) * (potential[index(i + 1, ny - 1, k)] + potential[index(i - 1, ny - 1, k)] +
                                                    potential[index(i, ny - 2, k)] + potential[index(i, ny - 1, k + 1)] +
                                                    potential[index(i, ny - 1, k - 1)]);
        }
    }

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            int idx = index(i, j, 0);
            if (!fixed_mask[idx])
                new_potential[idx] = (1.0 / 5.0) * (potential[index(i + 1, j, 0)] + potential[index(i - 1, j, 0)] +
                                                   potential[index(i, j + 1, 0)] + potential[index(i, j - 1, 0)] +
                                                   potential[index(i, j, 1)]);

            int idx1 = index(i, j, nz - 1);
            if (!fixed_mask[idx1])
                new_potential[idx1] = (1.0 / 5.0) * (potential[index(i + 1, j, nz - 1)] + potential[index(i - 1, j, nz - 1)] +
                                                    potential[index(i, j + 1, nz - 1)] + potential[index(i, j - 1, nz - 1)] +
                                                    potential[index(i, j, nz - 2)]);
        }
    }

    // Edge points (4 neighbors)
    for (int k = 1; k < nz - 1; ++k)
    {
        int idx = index(0, 0, k);
        if (!fixed_mask[idx])
            new_potential[idx] = (1.0 / 4.0) * (potential[index(1, 0, k)] + potential[index(0, 1, k)] +
                                               potential[index(0, 0, k + 1)] + potential[index(0, 0, k - 1)]);

        int idx1 = index(0, ny - 1, k);
        if (!fixed_mask[idx1])
            new_potential[idx1] = (1.0 / 4.0) * (potential[index(1, ny - 1, k)] + potential[index(0, ny - 2, k)] +
                                                potential[index(0, ny - 1, k + 1)] + potential[index(0, ny - 1, k - 1)]);

        int idx2 = index(nx - 1, 0, k);
        if (!fixed_mask[idx2])
            new_potential[idx2] = (1.0 / 4.0) * (potential[index(nx - 2, 0, k)] + potential[index(nx - 1, 1, k)] +
                                                potential[index(nx - 1, 0, k + 1)] + potential[index(nx - 1, 0, k - 1)]);

        int idx3 = index(nx - 1, ny - 1, k);
        if (!fixed_mask[idx3])
            new_potential[idx3] = (1.0 / 4.0) * (potential[index(nx - 2, ny - 1, k)] + potential[index(nx - 1, ny - 2, k)] +
                                                potential[index(nx - 1, ny - 1, k + 1)] + potential[index(nx - 1, ny - 1, k - 1)]);
    }

    for (int j = 1; j < ny - 1; ++j)
    {
        int idx = index(0, j, 0);
        if (!fixed_mask[idx])
            new_potential[idx] = (1.0 / 4.0) * (potential[index(1, j, 0)] + potential[index(0, j + 1, 0)] +
                                               potential[index(0, j - 1, 0)] + potential[index(0, j, 1)]);

        int idx1 = index(nx - 1, j, 0);
        if (!fixed_mask[idx1])
            new_potential[idx1] = (1.0 / 4.0) * (potential[index(nx - 2, j, 0)] + potential[index(nx - 1, j + 1, 0)] +
                                                potential[index(nx - 1, j - 1, 0)] + potential[index(nx - 1, j, 1)]);

        int idx2 = index(0, j, nz - 1);
        if (!fixed_mask[idx2])
            new_potential[idx2] = (1.0 / 4.0) * (potential[index(1, j, nz - 1)] + potential[index(0, j + 1, nz - 1)] +
                                                potential[index(0, j - 1, nz - 1)] + potential[index(0, j, nz - 2)]);

        int idx3 = index(nx - 1, j, nz - 1);
        if (!fixed_mask[idx3])
            new_potential[idx3] = (1.0 / 4.0) * (potential[index(nx - 2, j, nz - 1)] + potential[index(nx - 1, j + 1, nz - 1)] +
                                                potential[index(nx - 1, j - 1, nz - 1)] + potential[index(nx - 1, j, nz - 2)]);
    }

    for (int i = 1; i < nx - 1; ++i)
    {
        int idx = index(i, 0, 0);
        if (!fixed_mask[idx])
            new_potential[idx] = (1.0 / 4.0) * (potential[index(i + 1, 0, 0)] + potential[index(i - 1, 0, 0)] +
                                               potential[index(i, 1, 0)] + potential[index(i, 0, 1)]);

        int idx1 = index(i, ny - 1, 0);
        if (!fixed_mask[idx1])
            new_potential[idx1] = (1.0 / 4.0) * (potential[index(i + 1, ny - 1, 0)] + potential[index(i - 1, ny - 1, 0)] +
                                                potential[index(i, ny - 2, 0)] + potential[index(i, ny - 1, 1)]);

        int idx2 = index(i, 0, nz - 1);
        if (!fixed_mask[idx2])
            new_potential[idx2] = (1.0 / 4.0) * (potential[index(i + 1, 0, nz - 1)] + potential[index(i - 1, 0, nz - 1)] +
                                                potential[index(i, 1, nz - 1)] + potential[index(i, 0, nz - 2)]);

        int idx3 = index(i, ny - 1, nz - 1);
        if (!fixed_mask[idx3])
            new_potential[idx3] = (1.0 / 4.0) * (potential[index(i + 1, ny - 1, nz - 1)] + potential[index(i - 1, ny - 1, nz - 1)] +
                                                potential[index(i, ny - 2, nz - 1)] + potential[index(i, ny - 1, nz - 2)]);
    }

    // Corner points (3 neighbors)
    if (!fixed_mask[index(0, 0, 0)])
        new_potential[index(0, 0, 0)] = (1.0 / 3.0) * (potential[index(1, 0, 0)] + potential[index(0, 1, 0)] + potential[index(0, 0, 1)]);

    if (!fixed_mask[index(nx - 1, 0, 0)])
        new_potential[index(nx - 1, 0, 0)] = (1.0 / 3.0) * (potential[index(nx - 2, 0, 0)] + potential[index(nx - 1, 1, 0)] + potential[index(nx - 1, 0, 1)]);

    if (!fixed_mask[index(0, ny - 1, 0)])
        new_potential[index(0, ny - 1, 0)] = (1.0 / 3.0) * (potential[index(1, ny - 1, 0)] + potential[index(0, ny - 2, 0)] + potential[index(0, ny - 1, 1)]);

    if (!fixed_mask[index(nx - 1, ny - 1, 0)])
        new_potential[index(nx - 1, ny - 1, 0)] = (1.0 / 3.0) * (potential[index(nx - 2, ny - 1, 0)] + potential[index(nx - 1, ny - 2, 0)] + potential[index(nx - 1, ny - 1, 1)]);

    if (!fixed_mask[index(0, 0, nz - 1)])
        new_potential[index(0, 0, nz - 1)] = (1.0 / 3.0) * (potential[index(1, 0, nz - 1)] + potential[index(0, 1, nz - 1)] + potential[index(0, 0, nz - 2)]);

    if (!fixed_mask[index(nx - 1, 0, nz - 1)])
        new_potential[index(nx - 1, 0, nz - 1)] = (1.0 / 3.0) * (potential[index(nx - 2, 0, nz - 1)] + potential[index(nx - 1, 1, nz - 1)] + potential[index(nx - 1, 0, nz - 2)]);

    if (!fixed_mask[index(0, ny - 1, nz - 1)])
        new_potential[index(0, ny - 1, nz - 1)] = (1.0 / 3.0) * (potential[index(1, ny - 1, nz - 1)] + potential[index(0, ny - 2, nz - 1)] + potential[index(0, ny - 1, nz - 2)]);

    if (!fixed_mask[index(nx - 1, ny - 1, nz - 1)])
        new_potential[index(nx - 1, ny - 1, nz - 1)] = (1.0 / 3.0) * (potential[index(nx - 2, ny - 1, nz - 1)] + potential[index(nx - 1, ny - 2, nz - 1)] + potential[index(nx - 1, ny - 1, nz - 2)]);

    // Final update
    potential = std::move(new_potential);
}

void SimulationBox3D::applyGaussSeidel()
{
    // Interior points
    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            for (int k = 1; k < nz - 1; ++k)
            {
                int idx = index(i, j, k);
                if (!fixed_mask[idx])
                {
                    potential[idx] = (1.0 / 6.0) * (potential[index(i + 1, j, k)] + potential[index(i - 1, j, k)] +
                                                   potential[index(i, j + 1, k)] + potential[index(i, j - 1, k)] +
                                                   potential[index(i, j, k + 1)] + potential[index(i, j, k - 1)]);
                }
            }
        }
    }

    // Face points (excluding edges and corners)
    for (int j = 1; j < ny - 1; ++j)
    {
        for (int k = 1; k < nz - 1; ++k)
        {
            int idx = index(0, j, k);
            if (!fixed_mask[idx])
                potential[idx] = (1.0 / 5.0) * (potential[index(1, j, k)] + potential[index(0, j + 1, k)] +
                                               potential[index(0, j - 1, k)] + potential[index(0, j, k + 1)] +
                                               potential[index(0, j, k - 1)]);

            int idx1 = index(nx - 1, j, k);
            if (!fixed_mask[idx1])
                potential[idx1] = (1.0 / 5.0) * (potential[index(nx - 2, j, k)] + potential[index(nx - 1, j + 1, k)] +
                                                potential[index(nx - 1, j - 1, k)] + potential[index(nx - 1, j, k + 1)] +
                                                potential[index(nx - 1, j, k - 1)]);
        }
    }

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int k = 1; k < nz - 1; ++k)
        {
            int idx = index(i, 0, k);
            if (!fixed_mask[idx])
                potential[idx] = (1.0 / 5.0) * (potential[index(i + 1, 0, k)] + potential[index(i - 1, 0, k)] +
                                               potential[index(i, 1, k)] + potential[index(i, 0, k + 1)] +
                                               potential[index(i, 0, k - 1)]);

            int idx1 = index(i, ny - 1, k);
            if (!fixed_mask[idx1])
                potential[idx1] = (1.0 / 5.0) * (potential[index(i + 1, ny - 1, k)] + potential[index(i - 1, ny - 1, k)] +
                                                potential[index(i, ny - 2, k)] + potential[index(i, ny - 1, k + 1)] +
                                                potential[index(i, ny - 1, k - 1)]);
        }
    }

    for (int i = 1; i < nx - 1; ++i)
    {
        for (int j = 1; j < ny - 1; ++j)
        {
            int idx = index(i, j, 0);
            if (!fixed_mask[idx])
                potential[idx] = (1.0 / 5.0) * (potential[index(i + 1, j, 0)] + potential[index(i - 1, j, 0)] +
                                               potential[index(i, j + 1, 0)] + potential[index(i, j - 1, 0)] +
                                               potential[index(i, j, 1)]);

            int idx1 = index(i, j, nz - 1);
            if (!fixed_mask[idx1])
                potential[idx1] = (1.0 / 5.0) * (potential[index(i + 1, j, nz - 1)] + potential[index(i - 1, j, nz - 1)] +
                                                potential[index(i, j + 1, nz - 1)] + potential[index(i, j - 1, nz - 1)] +
                                                potential[index(i, j, nz - 2)]);
        }
    }

    // Edge points (4 neighbors)
    for (int k = 1; k < nz - 1; ++k)
    {
        // x edges at y=0 and y=ny-1
        int idx = index(0, 0, k);
        if (!fixed_mask[idx])
            potential[idx] = (1.0 / 4.0) * (potential[index(1, 0, k)] + potential[index(0, 1, k)] +
                                           potential[index(0, 0, k + 1)] + potential[index(0, 0, k - 1)]);

        int idx1 = index(0, ny - 1, k);
        if (!fixed_mask[idx1])
            potential[idx1] = (1.0 / 4.0) * (potential[index(1, ny - 1, k)] + potential[index(0, ny - 2, k)] +
                                            potential[index(0, ny - 1, k + 1)] + potential[index(0, ny - 1, k - 1)]);

        int idx2 = index(nx - 1, 0, k);
        if (!fixed_mask[idx2])
            potential[idx2] = (1.0 / 4.0) * (potential[index(nx - 2, 0, k)] + potential[index(nx - 1, 1, k)] +
                                            potential[index(nx - 1, 0, k + 1)] + potential[index(nx - 1, 0, k - 1)]);

        int idx3 = index(nx - 1, ny - 1, k);
        if (!fixed_mask[idx3])
            potential[idx3] = (1.0 / 4.0) * (potential[index(nx - 2, ny - 1, k)] + potential[index(nx - 1, ny - 2, k)] +
                                            potential[index(nx - 1, ny - 1, k + 1)] + potential[index(nx - 1, ny - 1, k - 1)]);
    }

    for (int j = 1; j < ny - 1; ++j)
    {
        int idx = index(0, j, 0);
        if (!fixed_mask[idx])
            potential[idx] = (1.0 / 4.0) * (potential[index(1, j, 0)] + potential[index(0, j + 1, 0)] +
                                           potential[index(0, j - 1, 0)] + potential[index(0, j, 1)]);

        int idx1 = index(nx - 1, j, 0);
        if (!fixed_mask[idx1])
            potential[idx1] = (1.0 / 4.0) * (potential[index(nx - 2, j, 0)] + potential[index(nx - 1, j + 1, 0)] +
                                            potential[index(nx - 1, j - 1, 0)] + potential[index(nx - 1, j, 1)]);

        int idx2 = index(0, j, nz - 1);
        if (!fixed_mask[idx2])
            potential[idx2] = (1.0 / 4.0) * (potential[index(1, j, nz - 1)] + potential[index(0, j + 1, nz - 1)] +
                                            potential[index(0, j - 1, nz - 1)] + potential[index(0, j, nz - 2)]);

        int idx3 = index(nx - 1, j, nz - 1);
        if (!fixed_mask[idx3])
            potential[idx3] = (1.0 / 4.0) * (potential[index(nx - 2, j, nz - 1)] + potential[index(nx - 1, j + 1, nz - 1)] +
                                            potential[index(nx - 1, j - 1, nz - 1)] + potential[index(nx - 1, j, nz - 2)]);
    }

    for (int i = 1; i < nx - 1; ++i)
    {
        int idx = index(i, 0, 0);
        if (!fixed_mask[idx])
            potential[idx] = (1.0 / 4.0) * (potential[index(i + 1, 0, 0)] + potential[index(i - 1, 0, 0)] +
                                           potential[index(i, 1, 0)] + potential[index(i, 0, 1)]);

        int idx1 = index(i, ny - 1, 0);
        if (!fixed_mask[idx1])
            potential[idx1] = (1.0 / 4.0) * (potential[index(i + 1, ny - 1, 0)] + potential[index(i - 1, ny - 1, 0)] +
                                            potential[index(i, ny - 2, 0)] + potential[index(i, ny - 1, 1)]);

        int idx2 = index(i, 0, nz - 1);
        if (!fixed_mask[idx2])
            potential[idx2] = (1.0 / 4.0) * (potential[index(i + 1, 0, nz - 1)] + potential[index(i - 1, 0, nz - 1)] +
                                            potential[index(i, 1, nz - 1)] + potential[index(i, 0, nz - 2)]);

        int idx3 = index(i, ny - 1, nz - 1);
        if (!fixed_mask[idx3])
            potential[idx3] = (1.0 / 4.0) * (potential[index(i + 1, ny - 1, nz - 1)] + potential[index(i - 1, ny - 1, nz - 1)] +
                                            potential[index(i, ny - 2, nz - 1)] + potential[index(i, ny - 1, nz - 2)]);
    }

    // Corners (3 neighbors)
    if (!fixed_mask[index(0, 0, 0)])
        potential[index(0, 0, 0)] = (1.0 / 3.0) * (potential[index(1, 0, 0)] + potential[index(0, 1, 0)] + potential[index(0, 0, 1)]);

    if (!fixed_mask[index(nx - 1, 0, 0)])
        potential[index(nx - 1, 0, 0)] = (1.0 / 3.0) * (potential[index(nx - 2, 0, 0)] + potential[index(nx - 1, 1, 0)] + potential[index(nx - 1, 0, 1)]);

    if (!fixed_mask[index(0, ny - 1, 0)])
        potential[index(0, ny - 1, 0)] = (1.0 / 3.0) * (potential[index(1, ny - 1, 0)] + potential[index(0, ny - 2, 0)] + potential[index(0, ny - 1, 1)]);

    if (!fixed_mask[index(nx - 1, ny - 1, 0)])
        potential[index(nx - 1, ny - 1, 0)] = (1.0 / 3.0) * (potential[index(nx - 2, ny - 1, 0)] + potential[index(nx - 1, ny - 2, 0)] + potential[index(nx - 1, ny - 1, 1)]);

    if (!fixed_mask[index(0, 0, nz - 1)])
        potential[index(0, 0, nz - 1)] = (1.0 / 3.0) * (potential[index(1, 0, nz - 1)] + potential[index(0, 1, nz - 1)] + potential[index(0, 0, nz - 2)]);

    if (!fixed_mask[index(nx - 1, 0, nz - 1)])
        potential[index(nx - 1, 0, nz - 1)] = (1.0 / 3.0) * (potential[index(nx - 2, 0, nz - 1)] + potential[index(nx - 1, 1, nz - 1)] + potential[index(nx - 1, 0, nz - 2)]);

    if (!fixed_mask[index(0, ny - 1, nz - 1)])
        potential[index(0, ny - 1, nz - 1)] = (1.0 / 3.0) * (potential[index(1, ny - 1, nz - 1)] + potential[index(0, ny - 2, nz - 1)] + potential[index(0, ny - 1, nz - 2)]);

    if (!fixed_mask[index(nx - 1, ny - 1, nz - 1)])
        potential[index(nx - 1, ny - 1, nz - 1)] = (1.0 / 3.0) * (potential[index(nx - 2, ny - 1, nz - 1)] + potential[index(nx - 1, ny - 2, nz - 1)] + potential[index(nx - 1, ny - 1, nz - 2)]);
}




/*
    CONVERSION OF CPP OBJECT TO PYTHON OBJECT
*/
// ---------- Convert C++ Value to PyObject ----------
template <typename T>
PyObject *cpp_value_to_pyobject(const T &value)
{
    if constexpr (std::is_same_v<T, int>)
        return PyLong_FromLong(value);
    else if constexpr (std::is_same_v<T, double>)
        return PyFloat_FromDouble(value);
    else if constexpr (std::is_same_v<T, float>)
        return PyFloat_FromDouble(static_cast<double>(value));
    else if constexpr (std::is_same_v<T, std::string>)
        return PyUnicode_FromString(value.c_str());
    else if constexpr (std::is_same_v<T, const char *>)
        return PyUnicode_FromString(value);
    else if constexpr (std::is_same_v<T, char *>)
        return PyUnicode_FromString(value);

    else if constexpr (std::is_same_v<T, char>)
        return PyUnicode_FromStringAndSize(&value, 1);
    else if constexpr (std::is_same_v<T, bool>)
        return PyBool_FromLong(value);
    else
    {
        static_assert(!std::is_same_v<T, T>, "Unsupported type for Python conversion.");
        return nullptr;
    }
}

/*
    CONVERSION OF SINGLE VARIABLE LIKE int, bool, char, double, string, ... etc
*/
// ---------- Inject Single Variable ----------
template <typename T>
void convert_and_inject_variable(const T &value, const char *name)
{
    PyObject *pyVal = cpp_value_to_pyobject(value);
    if (!pyVal)
        throw std::runtime_error("Failed to convert C++ value to Python object.");

    PyObject *mainModule = PyImport_AddModule("__main__");
    PyObject *mainDict = PyModule_GetDict(mainModule);
    PyDict_SetItemString(mainDict, name, pyVal);
    Py_DECREF(pyVal);
}

// ---------- Tuple Expansion ----------
template <typename TupleVars, typename TupleNames, std::size_t... Is>
void convert_and_inject_variables_impl(const TupleVars &vars, const TupleNames &names, std::index_sequence<Is...>)
{
    (convert_and_inject_variable(std::get<Is>(vars), std::get<Is>(names)), ...);
}

template <typename... Vars, typename... Names>
void convert_and_inject_variables(const std::tuple<Vars...> &vars, const std::tuple<Names...> &names)
{
    static_assert(sizeof...(Vars) == sizeof...(Names), "Mismatched number of variables and names.");
    convert_and_inject_variables_impl(vars, names, std::index_sequence_for<Vars...>{});
}

/*
    CONVERSION OF 1D VECTORS LIKE
    {1,2,6,4} -> [1,2,6,4]
    {1.5,-5.6,3e10} -> [1.5,-5.6,3e10]
    {TRUE, FALSE, FALSE,... } ->  [TRUE, FALSE, FALSE,... ]
*/
// ---------- Convert & Inject Vector ----------
template <typename T>
void convert_and_inject_vector(const std::vector<T> &vec, const char *name)
{
    PyObject *pyList = PyList_New(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
    {
        PyObject *item = cpp_value_to_pyobject(vec[i]);
        if (!item)
        {
            Py_DECREF(pyList);
            throw std::runtime_error("Failed to convert C++ type to Python object.");
        }
        PyList_SetItem(pyList, i, item); // Steals reference
    }

    PyObject *mainModule = PyImport_AddModule("__main__");
    PyObject *mainDict = PyModule_GetDict(mainModule);
    PyDict_SetItemString(mainDict, name, pyList);
    Py_DECREF(pyList); // PyDict_SetItemString does not steal
}

// ---------- Tuple Expansion ----------
template <typename TupleVecs, typename TupleNames, std::size_t... Is>
void convert_and_inject_vectors_impl(const TupleVecs &vecs, const TupleNames &names, std::index_sequence<Is...>)
{
    (convert_and_inject_vector(std::get<Is>(vecs), std::get<Is>(names)), ...);
}

template <typename... Vectors, typename... Names>
void convert_and_inject_vectors(const std::tuple<Vectors...> &vecs, const std::tuple<Names...> &names)
{
    static_assert(sizeof...(Vectors) == sizeof...(Names), "Number of vectors and names must match.");
    convert_and_inject_vectors_impl(vecs, names, std::index_sequence_for<Vectors...>{});
}

/*
    CONVERSION OF MULTIDIMENSIONAL VECTORS LIKE
    {{{1.1, 1.2}, {1.3, 1.4}}, {{2.1, 2.2}, {2.3, 2.4}}} -> [[[1.1 1.2],[1.3 1.4]],   [[2.1 2.2], [2.3 2.4]]]
    {{{true, false}, {false, true}}, {{true, true}, {false, false}}} ->  [[[True, False], [False, True]], [[True, True], [False, False]]]
*/
// Type trait to check if T is a std::vector
template <typename T>
struct is_std_vector : std::false_type
{
};

template <typename T, typename Alloc>
struct is_std_vector<std::vector<T, Alloc>> : std::true_type
{
};

template <typename T>
PyObject *cpp_tensorvalue_to_pyobject(const T &value)
{
    if constexpr (is_std_vector<T>::value)
    {
        PyObject *list = PyList_New(value.size());
        for (size_t i = 0; i < value.size(); ++i)
        {
            PyObject *item = cpp_value_to_pyobject(value[i]); // Recursive call
            PyList_SetItem(list, i, item);
        }
        return list;
    }
    else
    {
        // Scalar conversion
        if constexpr (std::is_same_v<T, int>)
            return PyLong_FromLong(value);
        else if constexpr (std::is_same_v<T, double>)
            return PyFloat_FromDouble(value);
        else if constexpr (std::is_same_v<T, bool>)
            return PyBool_FromLong(value);
        else if constexpr (std::is_same_v<T, const char *>)
            return PyUnicode_FromString(value);
        else if constexpr (std::is_same_v<T, std::string>)
            return PyUnicode_FromString(value.c_str());
        else
            static_assert(!sizeof(T), "Unsupported type for Python conversion");
    }
}

// Recursive converter
template <typename T>
PyObject *cpp_tensor_to_pyobject(const T &value)
{
    if constexpr (is_std_vector<T>::value)
    {
        PyObject *list = PyList_New(value.size());
        for (size_t i = 0; i < value.size(); ++i)
        {
            PyObject *item = cpp_tensor_to_pyobject(value[i]);
            if (!item)
            {
                Py_DECREF(list);
                throw std::runtime_error("Failed to convert nested vector.");
            }
            PyList_SetItem(list, i, item); // Steals reference
        }
        return list;
    }
    else
    {
        return cpp_tensorvalue_to_pyobject(value); // base case
    }
}

// Inject into Python
template <typename T>
void convert_and_inject_tensor(const T &tensor, const char *name)
{
    PyObject *pyTensor = cpp_tensor_to_pyobject(tensor);
    if (!pyTensor)
        throw std::runtime_error("Conversion failed.");

    PyObject *mainModule = PyImport_AddModule("__main__");
    PyObject *mainDict = PyModule_GetDict(mainModule);
    PyDict_SetItemString(mainDict, name, pyTensor);
    Py_DECREF(pyTensor);
}

// Tuple support
template <typename TupleTensors, typename TupleNames, std::size_t... Is>
void convert_and_inject_tensors_impl(const TupleTensors &tensors, const TupleNames &names, std::index_sequence<Is...>)
{
    (convert_and_inject_tensor(std::get<Is>(tensors), std::get<Is>(names)), ...);
}

template <typename... Tensors, typename... Names>
void convert_and_inject_tensors(const std::tuple<Tensors...> &tensors, const std::tuple<Names...> &names)
{
    static_assert(sizeof...(Tensors) == sizeof...(Names), "Mismatch in tensors and names.");
    convert_and_inject_tensors_impl(tensors, names, std::index_sequence_for<Tensors...>{});
}

/*
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

                                        MAIN FUNCTION

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
*/



int main() {
    // Define grid size and physical size
    int nx = 200, ny = 100, nz = 100;
    double lx = 2.0, ly = 1.0, lz = 1.0;

    SimulationBox3D box(nx, ny, nz, lx, ly, lz);

    // Add a sphere in the center at 1V
    // box.addSphere(0.5, 0.5, 0.5, 0.1, 1.0);

    // Add a box at a lower corner at -1V
    // box.addBox(0.1, 0.1, 0.1, 0.2, 0.2, 0.2, -1.0);

    box.addPlane(0,0,1,-0.3,0.01,2000);
    box.addPlane(0,0,1,-0.7,0.01,-2000);

    box.addCylinder(0,0.2,0.5,0.1,2,'x',2000);
    box.addCylinder(0,0.8,0.5,0.1,2,'x',-2000);



    // Solve the Laplace equation
    box.solve(5000, 5e-1, "gauss-seidel");

    // Check potential at center
    int i = nx / 2, j = ny / 2, k = nz / 2;
    double center_potential = box.getPotential()[box.index(i, j, k)];

    std::cout << "Potential at center: " << center_potential << " V" << std::endl;

    Py_Initialize();
    PyRun_SimpleString("from mayavi import mlab\nimport numpy as np");

    convert_and_inject_vectors(std::make_tuple(box.potential),
                                   std::make_tuple("box_potential"));
 
    convert_and_inject_variables(std::make_tuple(nx,ny,nz,lx,ly,lz),
                                    std::make_tuple("nx","ny","nz","lx","ly","lz"));

    const char *pyCode1 = R"(
box_potential = np.array(box_potential)
x, y, z = np.mgrid[0:lx:nx*1j,0:ly:ny*1j,0:lz:nz*1j]
contours=100
opacity=0.08
cmap="jet"
print(x.shape)
box_potential=box_potential.reshape(x.shape)
print(box_potential.shape)

mlab.figure(bgcolor=(1, 1, 1), size=(800, 600))
mlab.contour3d(x, y, z, box_potential, contours=contours, opacity=opacity, colormap=cmap)
axes = mlab.axes(xlabel='X', ylabel='Y', zlabel='Z',color=(1.0,0.0,0.0))
mlab.colorbar(title="potential", orientation='vertical')
axes.title_text_property.color = (1.0, 0.0, 0.0)  # red title text (if titles used)
axes.label_text_property.color = (1.0, 0.0, 0.0)  # blue label text
mlab.title("3D Isosurface of Potential")
mlab.show()

        )";

    PyRun_SimpleString(pyCode1);


    return 0;
}