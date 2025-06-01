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

double Bx = 5;
double By = 0;
double Bz = 0;

double cm = 1e-2;
double mm = 1e-3;
double qe = 1.60217663e-19;
double me = 9.10938e-31;
double mH = 1.67262192e-27;

double kev_to_joule = 1.660217663e-16;

/*
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

                                3D SIMULATION BOX OBJECT FOR MAIN GEOMETRY

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

class SimulationBox3D
{
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

    void solve(int max_iter = 1000, double tol = 1e-4, const std::string &method = "jacobi");

    // Add geometry will be added next (e.g., addBox, addSphere)

    const std::vector<double> &getPotential() const { return potential; }

    std::vector<double> potential;
    std::vector<bool> fixed_mask;

    int index(int i, int j, int k) const
    {
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

void SimulationBox3D::solve(int max_iter, double tol, const std::string &method)
{
    std::vector<double> V_old = potential;
    std::vector<double> V_new = potential;

    for (int iter = 0; iter < max_iter; ++iter)
    {
        V_old = potential;

        if (method == "jacobi")
        {
            applyJacobi();
        }
        else if (method == "gauss-seidel")
        {
            applyGaussSeidel();
        }
        else
        {
            throw std::runtime_error("Unknown method");
        }

        // Check for convergence
        double max_diff = 0.0;
        for (size_t i = 0; i < potential.size(); ++i)
        {
            max_diff = std::max(max_diff, std::abs(potential[i] - V_old[i]));
        }

        std::cout << "Iteration " << iter << ", Max Diff = " << max_diff << std::endl;

        if (max_diff < tol)
        {
            std::cout << "Converged at iteration " << iter << std::endl;
            break;
        }
    }
}

void SimulationBox3D::addSphere(double cx, double cy, double cz, double radius, double potential_value)
{
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        for (int j = 0; j < ny; ++j)
        {
            double y = j * dy;
            for (int k = 0; k < nz; ++k)
            {
                double z = k * dz;
                double dist = std::sqrt((x - cx) * (x - cx) +
                                        (y - cy) * (y - cy) +
                                        (z - cz) * (z - cz));
                if (dist <= radius)
                {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}

void SimulationBox3D::addBox(double x0, double y0, double z0,
                             double x1, double y1, double z1, double potential_value)
{
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        if (x < x0 || x > x1)
            continue;
        for (int j = 0; j < ny; ++j)
        {
            double y = j * dy;
            if (y < y0 || y > y1)
                continue;
            for (int k = 0; k < nz; ++k)
            {
                double z = k * dz;
                if (z < z0 || z > z1)
                    continue;

                int idx = index(i, j, k);
                potential[idx] = potential_value;
                fixed_mask[idx] = true;
            }
        }
    }
}

void SimulationBox3D::addCylinder(double cx, double cy, double cz, double radius, double height, char axis, double potential_value)
{
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        for (int j = 0; j < ny; ++j)
        {
            double y = j * dy;
            for (int k = 0; k < nz; ++k)
            {
                double z = k * dz;

                bool inside = false;
                if (axis == 'z')
                {
                    double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                    inside = (r2 <= radius * radius && z >= cz && z <= cz + height);
                }
                else if (axis == 'x')
                {
                    double r2 = (y - cy) * (y - cy) + (z - cz) * (z - cz);
                    inside = (r2 <= radius * radius && x >= cx && x <= cx + height);
                }
                else if (axis == 'y')
                {
                    double r2 = (x - cx) * (x - cx) + (z - cz) * (z - cz);
                    inside = (r2 <= radius * radius && y >= cy && y <= cy + height);
                }

                if (inside)
                {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}

void SimulationBox3D::addHollowPipe(double cx, double cy, double cz, double radius, double thickness, double height, char axis, double potential_value)
{
    double r_outer2 = (radius + thickness / 2.0) * (radius + thickness / 2.0);
    double r_inner2 = (radius - thickness / 2.0) * (radius - thickness / 2.0);

    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        for (int j = 0; j < ny; ++j)
        {
            double y = j * dy;
            for (int k = 0; k < nz; ++k)
            {
                double z = k * dz;

                bool inside = false;
                if (axis == 'z')
                {
                    double r2 = (x - cx) * (x - cx) + (y - cy) * (y - cy);
                    inside = (r2 >= r_inner2 && r2 <= r_outer2 && z >= cz && z <= cz + height);
                }
                else if (axis == 'x')
                {
                    double r2 = (y - cy) * (y - cy) + (z - cz) * (z - cz);
                    inside = (r2 >= r_inner2 && r2 <= r_outer2 && x >= cx && x <= cx + height);
                }
                else if (axis == 'y')
                {
                    double r2 = (x - cx) * (x - cx) + (z - cz) * (z - cz);
                    inside = (r2 >= r_inner2 && r2 <= r_outer2 && y >= cy && y <= cy + height);
                }

                if (inside)
                {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}

void SimulationBox3D::addEllipsoid(double cx, double cy, double cz, double rx, double ry, double rz, double potential_value)
{
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        for (int j = 0; j < ny; ++j)
        {
            double y = j * dy;
            for (int k = 0; k < nz; ++k)
            {
                double z = k * dz;
                double value = ((x - cx) / rx) * ((x - cx) / rx) + ((y - cy) / ry) * ((y - cy) / ry) + ((z - cz) / rz) * ((z - cz) / rz);
                if (value <= 1.0)
                {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}

void SimulationBox3D::addHyperboloid(double cx, double cy, double cz, double a, double b, double c, double waist, char axis, double potential_value)
{
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        for (int j = 0; j < ny; ++j)
        {
            double y = j * dy;
            for (int k = 0; k < nz; ++k)
            {
                double z = k * dz;
                double val = 0.0;
                if (axis == 'x')
                {
                    val = -((x - cx) / a) * ((x - cx) / a) + ((y - cy) / b) * ((y - cy) / b) + ((z - cz) / c) * ((z - cz) / c);
                }
                else if (axis == 'y')
                {
                    val = ((x - cx) / a) * ((x - cx) / a) - ((y - cy) / b) * ((y - cy) / b) + ((z - cz) / c) * ((z - cz) / c);
                }
                else if (axis == 'z')
                {
                    val = ((x - cx) / a) * ((x - cx) / a) + ((y - cy) / b) * ((y - cy) / b) - ((z - cz) / c) * ((z - cz) / c);
                }

                if (val <= waist * waist)
                {
                    int idx = index(i, j, k);
                    potential[idx] = potential_value;
                    fixed_mask[idx] = true;
                }
            }
        }
    }
}

void SimulationBox3D::addPlane(double A, double B, double C, double D, double thickness, double potential_value)
{
    double norm = std::sqrt(A * A + B * B + C * C);
    for (int i = 0; i < nx; ++i)
    {
        double x = i * dx;
        for (int j = 0; j < ny; ++j)
        {
            double y = j * dy;
            for (int k = 0; k < nz; ++k)
            {
                double z = k * dz;
                double dist = (A * x + B * y + C * z + D) / norm;
                if (std::abs(dist) <= thickness / 2.0)
                {
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
    potential = std::move(new_potential);
}

void SimulationBox3D::applyGaussSeidel()
{
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
}

/*
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

                                DEFINE PARTICLE CLASS TO MOVE IN FIELD

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

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

void propagator(Particle &p,
                const SimulationBox3D &box,
                double t_max = 0.0,
                double dt = 0.001) // constant B field
{
    double t = 0.0;
    int steps = static_cast<int>(t_max / dt);
    double qmdt2 = (p.q / p.m) * (dt / 2.0);

    for (int step = 0; step < steps; ++step)
    {
        // Boundary check
        if (p.x < 0 || p.x >= box.lx ||
            p.y < 0 || p.y >= box.ly ||
            p.z < 0 || p.z >= box.lz)
        {
            std::cout << "Particle left the domain at step " << step << std::endl;
            break;
        }

        // Convert position to indices
        int i = static_cast<int>(p.x / box.dx);
        int j = static_cast<int>(p.y / box.dy);
        int k = static_cast<int>(p.z / box.dz);

        // Check if inside grid bounds
        if (i < 0 || i >= box.nx ||
            j < 0 || j >= box.ny ||
            k < 0 || k >= box.nz)
        {
            std::cout << "Particle reached grid edge at step " << step << std::endl;
            break;
        }

        // Check if on electrode
        if (box.fixed_mask[box.index(i, j, k)])
        {
            std::cout << "Particle hit an electrode at step " << step << std::endl;
            break;
        }

        // Interpolate E field (simple nearest-neighbor for now)
        double Ex = 0.0, Ey = 0.0, Ez = 0.0;

        // Use finite differences to approximate E = -∇φ
        if (i > 0 && i < box.nx - 1)
            Ex = -(box.potential[box.index(i + 1, j, k)] - box.potential[box.index(i - 1, j, k)]) / (2.0 * box.dx);
        if (j > 0 && j < box.ny - 1)
            Ey = -(box.potential[box.index(i, j + 1, k)] - box.potential[box.index(i, j - 1, k)]) / (2.0 * box.dy);
        if (k > 0 && k < box.nz - 1)
            Ez = -(box.potential[box.index(i, j, k + 1)] - box.potential[box.index(i, j, k - 1)]) / (2.0 * box.dz);

        // Half-step velocity
        double vx_minus = p.vx + qmdt2 * Ex;
        double vy_minus = p.vy + qmdt2 * Ey;
        double vz_minus = p.vz + qmdt2 * Ez;

        // Magnetic rotation
        double tx = qmdt2 * Bx;
        double ty = qmdt2 * By;
        double tz = qmdt2 * Bz;

        double t2 = tx * tx + ty * ty + tz * tz;
        double sx = 2 * tx / (1 + t2);
        double sy = 2 * ty / (1 + t2);
        double sz = 2 * tz / (1 + t2);

        double vpx = vx_minus + (vy_minus * tz - vz_minus * ty);
        double vpy = vy_minus + (vz_minus * tx - vx_minus * tz);
        double vpz = vz_minus + (vx_minus * ty - vy_minus * tx);

        double vx_plus = vx_minus + (vpy * sz - vpz * sy);
        double vy_plus = vy_minus + (vpz * sx - vpx * sz);
        double vz_plus = vz_minus + (vpx * sy - vpy * sx);

        // Final update
        p.vx = vx_plus + qmdt2 * Ex;
        p.vy = vy_plus + qmdt2 * Ey;
        p.vz = vz_plus + qmdt2 * Ez;

        // Position update
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;

        // Save trajectory
        p.posx.push_back(p.x);
        p.posy.push_back(p.y);
        p.posz.push_back(p.z);
    }

    // Final stats
    p.v = std::sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
    p.energy = 0.5 * p.m * p.v * p.v;
}

/*
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

                                CONVERSION OF CPP OBJECT TO PYTHON OBJECT

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
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

int main()
{
    // Define grid size and physical size
    int nx = 200, ny = 200, nz = 200;
    double lx = 2.0 * cm, ly = 1.0 * cm, lz = 1.0 * cm;

    SimulationBox3D box(nx, ny, nz, lx, ly, lz);

    // Add a sphere in the center at 1V
    // box.addSphere(0.5, 0.5, 0.5, 0.1, 1.0);

    // Add a box at a lower corner at -1V
    // box.addBox(0.1, 0.1, 0.1, 0.2, 0.2, 0.2, -1.0);

    box.addPlane(0, 0, 1, -0.3 * cm, 0.03 * cm, 2000);
    box.addPlane(0, 0, 1, -0.7 * cm, 0.03 * cm, 2000);

    box.addCylinder(0 * cm, 0.2 * cm, 0.5 * cm, 0.05 * cm, 2 * cm, 'x', 2000);
    box.addCylinder(0 * cm, 0.8 * cm, 0.5 * cm, 0.05 * cm, 2 * cm, 'x', -2000);

    // Solve the Laplace equation
    box.solve(5000, 2, "gauss-seidel");

    // Check potential at center
    int i = nx / 2, j = ny / 2, k = nz / 2;
    double center_potential = box.getPotential()[box.index(i, j, k)];

    std::cout << "Potential at center: " << center_potential << " V" << std::endl;

    double energy = 0.1 * kev_to_joule;
    double vx = sqrt(2 * energy / mH);

    Particle p(0.1 * cm, 0.5 * cm, 0.5 * cm, vx, -vx/2, -vx/2, qe, mH); // x, y, z, vx, vy, vz, q, m
    int res_t = 20000;
    double t_max = 4.0 * cm / p.v;
    double dt = t_max / res_t;
    propagator(p, box, t_max, dt);

    Py_Initialize();
    PyRun_SimpleString("from mayavi import mlab\nimport numpy as np");

    convert_and_inject_vectors(std::make_tuple(box.potential),
                               std::make_tuple("box_potential"));

    convert_and_inject_variables(std::make_tuple(nx, ny, nz, lx, ly, lz),
                                 std::make_tuple("nx", "ny", "nz", "lx", "ly", "lz"));

    convert_and_inject_vectors(std::make_tuple(p.posx, p.posy, p.posz),
                               std::make_tuple("p_posx", "p_posy", "p_posz"));

    const char *pyCode1 = R"(
box_potential = np.array(box_potential)
p_posx=np.array(p_posx)
p_posy=np.array(p_posy)
p_posz=np.array(p_posz)

print(p_posx)
print(p_posy)
print(p_posz)

x, y, z = np.mgrid[0:lx:nx*1j,0:ly:ny*1j,0:lz:nz*1j]
contours=50
opacity=0.08
cmap="jet"
print(x.shape)
box_potential=box_potential.reshape(x.shape)
print(box_potential.shape)

mlab.figure(bgcolor=(1, 1, 1), size=(800, 600))
mlab.contour3d(x, y, z, box_potential, contours=contours, opacity=opacity, colormap=cmap)
mlab.colorbar(title="potential", orientation='vertical')
axes = mlab.axes(xlabel='X', ylabel='Y', zlabel='Z',color=(1.0,0.0,0.0))
axes.title_text_property.color = (1.0, 0.0, 0.0)  # red title text (if titles used)
axes.label_text_property.color = (1.0, 0.0, 0.0)  # blue label text
mlab.title("3D Isosurface of Potential")
mlab.plot3d(p_posx,p_posy,p_posz, tube_radius=0.005*1e-2, color=(0,0,0), tube_sides=12)
mlab.show()
)";

    PyRun_SimpleString(pyCode1);

    return 0;
}