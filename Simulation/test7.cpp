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


// Define a particle class
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


// propagator function for the particle to move it in the EM field 
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


// function for printing a 1-d vector
void print_1d_vector(const std::vector<double> &v)
{
    for (const auto &elem : v)
    {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

// function to concatenate double vector like np.concatenate. y=concatenate(x,y,z)
// will  give [x1,x2,x3,.... , y1,y2,y3, .... , z1,z2,z3,....] 
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


// make a linspace for double value like np.linspace
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


// solve a 2x2 matrix equation (helper functon)
std::pair<double, double> solve2x2(double a11, double a12, double a21, double a22, double b1, double b2)
{
    double det = a11 * a22 - a12 * a21;
    if (std::abs(det) < 1e-12)
        throw std::runtime_error("Singular matrix");

    double alpha = (b1 * a22 - b2 * a12) / det;
    double beta = (a11 * b2 - a21 * b1) / det;

    return {alpha, beta};
}


// get a exponential distribution for particle with given steepness
// more particle in low energy and less in high
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


/*
    CONVERSION OF CPP OBJECT TO PYTHON OBJECT
*/
// ---------- Convert C++ Value to PyObject ----------
template <typename T>
PyObject* cpp_value_to_pyobject(const T& value)
{
    if constexpr (std::is_same_v<T, int>)
        return PyLong_FromLong(value);
    else if constexpr (std::is_same_v<T, double>)
        return PyFloat_FromDouble(value);
    else if constexpr (std::is_same_v<T, float>)
        return PyFloat_FromDouble(static_cast<double>(value));
    else if constexpr (std::is_same_v<T, std::string>)
        return PyUnicode_FromString(value.c_str());
    else if constexpr (std::is_same_v<T, const char*>)
        return PyUnicode_FromString(value);
    else if constexpr (std::is_same_v<T, char*>)
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
void convert_and_inject_variable(const T& value, const char* name)
{
    PyObject* pyVal = cpp_value_to_pyobject(value);
    if (!pyVal)
        throw std::runtime_error("Failed to convert C++ value to Python object.");

    PyObject* mainModule = PyImport_AddModule("__main__");
    PyObject* mainDict = PyModule_GetDict(mainModule);
    PyDict_SetItemString(mainDict, name, pyVal);
    Py_DECREF(pyVal);
}

// ---------- Tuple Expansion ----------
template <typename TupleVars, typename TupleNames, std::size_t... Is>
void convert_and_inject_variables_impl(const TupleVars& vars, const TupleNames& names, std::index_sequence<Is...>)
{
    (convert_and_inject_variable(std::get<Is>(vars), std::get<Is>(names)), ...);
}

template <typename... Vars, typename... Names>
void convert_and_inject_variables(const std::tuple<Vars...>& vars, const std::tuple<Names...>& names)
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
void convert_and_inject_vector(const std::vector<T>& vec, const char* name)
{
    PyObject* pyList = PyList_New(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
    {
        PyObject* item = cpp_value_to_pyobject(vec[i]);
        if (!item)
        {
            Py_DECREF(pyList);
            throw std::runtime_error("Failed to convert C++ type to Python object.");
        }
        PyList_SetItem(pyList, i, item); // Steals reference
    }

    PyObject* mainModule = PyImport_AddModule("__main__");
    PyObject* mainDict = PyModule_GetDict(mainModule);
    PyDict_SetItemString(mainDict, name, pyList);
    Py_DECREF(pyList);  // PyDict_SetItemString does not steal
}



// ---------- Tuple Expansion ----------
template <typename TupleVecs, typename TupleNames, std::size_t... Is>
void convert_and_inject_vectors_impl(const TupleVecs& vecs, const TupleNames& names, std::index_sequence<Is...>)
{
    (convert_and_inject_vector(std::get<Is>(vecs), std::get<Is>(names)), ...);
}



template <typename... Vectors, typename... Names>
void convert_and_inject_vectors(const std::tuple<Vectors...>& vecs, const std::tuple<Names...>& names)
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
PyObject* cpp_tensorvalue_to_pyobject(const T& value) {
    if constexpr (is_std_vector<T>::value) {
        PyObject* list = PyList_New(value.size());
        for (size_t i = 0; i < value.size(); ++i) {
            PyObject* item = cpp_value_to_pyobject(value[i]);  // Recursive call
            PyList_SetItem(list, i, item);
        }
        return list;
    } else {
        // Scalar conversion
        if constexpr (std::is_same_v<T, int>)
            return PyLong_FromLong(value);
        else if constexpr (std::is_same_v<T, double>)
            return PyFloat_FromDouble(value);
        else if constexpr (std::is_same_v<T, bool>)
            return PyBool_FromLong(value);
        else if constexpr (std::is_same_v<T, const char*>)
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

    // C++ variables of different types
    int i = 42;
    double d = 3.14159;
    float f = 2.71f;
    std::string s = "v_string";
    const char* cstr = "v_const_char*";
    char* cs ="v_char*";   // this will give warning... use const char* or string instead. char* is bad for memory
    char ch = 'X';
    bool flag = true;


    std::vector<int> vi = {1, 2, 3, 4, 5, 6};
    std::vector<float> vf = {1.2f, 2.3f, -3.56f, 4.451f, -50.365f, 6.0f};
    std::vector<double> vD = {10.13, -2.131, 3.23e12, -4.1e-10, 5.23, 60465.465};
    std::vector<char> vc = {'a', 'n', 'h', 'r'};
    std::vector<std::string> vs = {"v_string", "6546", "sjdh", "jd", "ajdhbajsjk", "5616as5a_$%$"};
    std::vector<const char*> vcstr = {"v_const_char*", "6546", "sjdh", "jd", "ajdhbajsjk", "5616as5a_$%$"};  
    std::vector<char*> vcs = {"v_char*", "6546", "sjdh", "jd", "ajdhbajsjk", "5616as5a_$%$"}; // this will give warning... use const char* or string instead. char* is bad for memory
    std::vector<bool> vb = {true, false, true, true};

    // Inject into Python
    convert_and_inject_variables(
        std::make_tuple(i, d, f, s, cstr, ch, flag, cs),
        std::make_tuple("i", "d", "f", "s", "cstr", "ch", "flag", "cs")
    );


    convert_and_inject_vectors(
        std::make_tuple(vi,vf,vD,vc,vs,vcs,vb,vcstr),
        std::make_tuple("vi", "vf", "vD", "vc", "vs", "vcs", "vb", "vcstr")
    );

    // Run a Python script that prints them
    const char* pyScript = R"(
print(f"i     (int)   = {i}")
print(f"d     (double)= {d}")
print(f"f     (float) = {f}")
print(f"s     (str)   = {s}")
print(f"cstr  (str)   = {cstr}")
print(f"char* (char*) = {cs}")
print(f"ch    (char)  = {ch}")
print(f"flag  (bool)  = {flag}")

print("Arrays\n___________")
print(f"vi: {vi}\nvf: {vf}\nvD: {vD}\nvc: {vc}\nvs: {vs}\nvcs: {vcs}\nvcstr: {vcstr}\nvb: {vb}")
)";

    PyRun_SimpleString(pyScript);

    std::vector<std::vector<std::vector<double>>> tensor3d = {
        {{1.1, 1.2}, {1.3, 1.4}},
        {{2.1, 2.2}, {2.3, 2.4}}};

    std::vector<std::vector<std::vector<bool>>> bool_tensor = {
        {{true, false}, {false, true}},
        {{true, true}, {false, false}}};

    convert_and_inject_tensors(
        std::make_tuple(tensor3d, bool_tensor),
        std::make_tuple("tensor3d", "bool_tensor"));

    PyRun_SimpleString(R"(
import numpy as np
print("3D Tensor:")
print(np.array(tensor3d))
print("\n3D Bool Tensor:")
print(bool_tensor)
)");

    Py_Finalize();
    return 0;

}
