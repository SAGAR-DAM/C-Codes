#include "C:\\Users\\mrsag\\AppData\\Local\\Programs\\Python\\Python311\\include\\Python.h"
#include <vector>
#include <tuple>
#include <iostream>

// Converts one vector and injects it into the __main__ module under a given name
void convert_and_inject_vector(const std::vector<double>& vec, const char* name) {
    PyObject* pyList = PyList_New(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        PyList_SetItem(pyList, i, PyFloat_FromDouble(vec[i]));
    }

    PyObject* mainModule = PyImport_AddModule("__main__");
    PyObject* mainDict = PyModule_GetDict(mainModule);
    PyDict_SetItemString(mainDict, name, pyList);

    Py_DECREF(pyList);
}

// Helper function to expand vector-name pairs by index
template <typename TupleVecs, typename TupleNames, std::size_t... Is>
void convert_and_inject_vectors_impl(const TupleVecs& vecs, const TupleNames& names, std::index_sequence<Is...>) {
    (convert_and_inject_vector(std::get<Is>(vecs), std::get<Is>(names)), ...);
}

// Main function: accepts tuples of vectors and names
template <typename... Vectors, typename... Names>
void convert_and_inject_vectors(const std::tuple<Vectors...>& vecs, const std::tuple<Names...>& names) {
    static_assert(sizeof...(Vectors) == sizeof...(Names), "Number of vectors and names must match.");
    convert_and_inject_vectors_impl(vecs, names, std::index_sequence_for<Vectors...>{});
}



int main() {
    std::vector<double> x = {0, 1, 2, 3, 4, 5};
    std::vector<double> y = {0, 1, 0, 1, 0, 1};
    std::vector<double> w = {5, 4, 3, 2, 1, 0};

    Py_Initialize();

    PyRun_SimpleString("import matplotlib.pyplot as plt\nimport numpy as np");

    // Use tuples for the vectors and names
    convert_and_inject_vectors(std::make_tuple(x, y, w),
                               std::make_tuple("x", "y", "w"));

    const char* pyCode = R"(
plt.figure()
plt.plot(x, y, label='x vs y')
plt.plot(x, w, label='x vs w')
plt.legend()
plt.grid(True)
plt.show()
)";
    PyRun_SimpleString(pyCode);

    Py_Finalize();
    return 0;
}
