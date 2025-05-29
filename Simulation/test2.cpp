#include "C:\\Users\\mrsag\\AppData\\Local\\Programs\\Python\\Python311\\include\\Python.h"
#include <iostream>

int main() {
    float a = 305.545f;

    Py_Initialize();

    // Create a Python float object from C++ float
    PyObject* pyMainModule = PyImport_AddModule("__main__");
    PyObject* pyGlobalDict = PyModule_GetDict(pyMainModule);
    PyDict_SetItemString(pyGlobalDict, "a", PyFloat_FromDouble(a));  // Inject C++ float 'a' into Python globals

    const char* pyCode1 = R"(
print("Value of a from C++:", a)

import matplotlib.pyplot as plt
plt.grid(True)
plt.plot([0], [0], "ko", markersize=5)
plt.gca().set_aspect('equal')
plt.show()
    )";

    PyRun_SimpleString(pyCode1);

    Py_Finalize();
    return 0;
}
