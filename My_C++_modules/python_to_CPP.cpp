#include <Python.h>
#include <iostream>

int main() {
    Py_Initialize();

    // Create a main Python namespace (like running a script)
    PyObject* main_module = PyImport_AddModule("__main__");
    PyObject* global_dict = PyModule_GetDict(main_module);

    // Set some variable from C++
    PyRun_SimpleString("x = 2\n"
                       "y = [1, 2, 3]\n"
                       "z = [[1, 2], [3, 4]]\n"
                       "s = 'hello'\n");

    // Get variable `x` from Python
    PyObject* pyX = PyDict_GetItemString(global_dict, "x");
    if (pyX && PyLong_Check(pyX)) {
        long x_value = PyLong_AsLong(pyX);
        std::cout << "x = " << x_value << std::endl;
    }

    // Get list `y` from Python
    PyObject* pyY = PyDict_GetItemString(global_dict, "y");
    if (pyY && PyList_Check(pyY)) {
        Py_ssize_t size = PyList_Size(pyY);
        std::cout << "y = [";
        for (Py_ssize_t i = 0; i < size; ++i) {
            PyObject* item = PyList_GetItem(pyY, i);
            long val = PyLong_AsLong(item);
            std::cout << val;
            if (i < size - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }

    // Get 2D list `z` from Python
    PyObject* pyZ = PyDict_GetItemString(global_dict, "z");
    if (pyZ && PyList_Check(pyZ)) {
        std::cout << "z = [\n";
        for (Py_ssize_t i = 0; i < PyList_Size(pyZ); ++i) {
            PyObject* inner = PyList_GetItem(pyZ, i);
            std::cout << "  [";
            for (Py_ssize_t j = 0; j < PyList_Size(inner); ++j) {
                PyObject* item = PyList_GetItem(inner, j);
                long val = PyLong_AsLong(item);
                std::cout << val;
                if (j < PyList_Size(inner) - 1) std::cout << ", ";
            }
            std::cout << "]\n";
        }
        std::cout << "]" << std::endl;
    }

    // Get string `s` from Python
    PyObject* pyS = PyDict_GetItemString(global_dict, "s");
    if (pyS && PyUnicode_Check(pyS)) {
        const char* str = PyUnicode_AsUTF8(pyS);
        std::cout << "s = " << str << std::endl;
    }

    Py_Finalize();
    return 0;
}
