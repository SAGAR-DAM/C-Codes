#include <Python.h>
#include <iostream>
#include <tuple>
#include <vector>
#include <string>
#include <type_traits>

// ---------- Convert C++ Value to PyObject ----------

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
    CONVERSION OF MULTIDIMENSIONAL VECTORS TO PYTHON ARRAY LIKE
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

// ---------- Main ----------
int main()
{
    Py_Initialize();

    // C++ variables of different types
    int i = 42;
    double d = 3.14159;
    float f = 2.71f;
    std::string s = "v_string";
    const char *cstr = "v_const_char*";
    char *cs = "v_char*"; // this will give warning... use const char* or string instead. char* is bad for memory
    char ch = 'X';
    bool flag = true;

    std::vector<int> vi = {1, 2, 3, 4, 5, 6};
    std::vector<float> vf = {1.2f, 2.3f, -3.56f, 4.451f, -50.365f, 6.0f};
    std::vector<double> vD = {10.13, -2.131, 3.23e12, -4.1e-10, 5.23, 60465.465};
    std::vector<char> vc = {'a', 'n', 'h', 'r'};
    std::vector<std::string> vs = {"v_string", "6546", "sjdh", "jd", "ajdhbajsjk", "5616as5a_$%$"};
    std::vector<const char *> vcstr = {"v_const_char*", "6546", "sjdh", "jd", "ajdhbajsjk", "5616as5a_$%$"};
    std::vector<char *> vcs = {"v_char*", "6546", "sjdh", "jd", "ajdhbajsjk", "5616as5a_$%$"}; // this will give warning... use const char* or string instead. char* is bad for memory
    std::vector<bool> vb = {true, false, true, true};

    // Inject into Python
    convert_and_inject_variables(
        std::make_tuple(i, d, f, s, cstr, ch, flag, cs),
        std::make_tuple("i", "d", "f", "s", "cstr", "ch", "flag", "cs"));

    convert_and_inject_vectors(
        std::make_tuple(vi, vf, vD, vc, vs, vcs, vb, vcstr),
        std::make_tuple("vi", "vf", "vD", "vc", "vs", "vcs", "vb", "vcstr"));

    // Run a Python script that prints them
    const char *pyScript = R"(
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
