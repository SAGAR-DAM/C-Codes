// AvgModule.cpp
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "C:\Users\mrsag\AppData\Local\Programs\Python\Python311\include\Python.h"
#include "C:\Users\mrsag\AppData\Local\Programs\Python\Python311\Lib\site-packages\numpy\core\include\numpy\arrayobject.h"

// Run with:
// g++ -std=c++11 AvgModule.cpp `
// -I"C:\Users\mrsag\AppData\Local\Programs\Python\Python311\include" `
// -I"C:\Users\mrsag\AppData\Local\Programs\Python\Python311\Lib\site-packages\numpy\core\include" `
// -L"C:\Users\mrsag\AppData\Local\Programs\Python\Python311\libs" `
// -lpython311 -shared -o AvgModule.pyd


// and if want to run with python 3.12: 
// g++ -std=c++11 AvgModule.cpp `
//   -I"C:\Users\mrsag\anaconda3\include" `
//   -I"C:\Users\mrsag\anaconda3\Lib\site-packages\numpy\core\include" `
//   -L"C:\Users\mrsag\anaconda3\libs" `
//   -lpython312 -shared -o AvgModule.pyd


static PyObject* cpp_average(PyObject* self, PyObject* args) {
    PyArrayObject* input_array;

    // Parse NumPy array from Python args
    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &input_array))
        return NULL;

    if (PyArray_TYPE(input_array) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "Only float64 arrays are supported.");
        return NULL;
    }

    double* data = (double*) PyArray_DATA(input_array);
    npy_intp size = PyArray_SIZE(input_array);

    double sum = 0.0;
    for (npy_intp i = 0; i < size; ++i) {
        sum += data[i];
    }

    double avg = sum / size;

    return Py_BuildValue("d", avg);
}

// Define methods
static PyMethodDef AvgMethods[] = {
    {"cpp_average", cpp_average, METH_VARARGS, "Compute average of a NumPy array"},
    {NULL, NULL, 0, NULL}
};

// Define module
static struct PyModuleDef AvgModule = {
    PyModuleDef_HEAD_INIT,
    "AvgModule",
    "Module for computing average using C++",
    -1,
    AvgMethods
};

// Init function
PyMODINIT_FUNC PyInit_AvgModule(void) {
    import_array();  // Initialize NumPy API
    return PyModule_Create(&AvgModule);
}
