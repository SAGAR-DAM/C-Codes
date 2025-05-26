#include "C:\Users\mrsag\AppData\Local\Programs\Python\Python311\Lib\site-packages\pybind11\include\pybind11\pybind11.h"
#include "C:\Users\mrsag\AppData\Local\Programs\Python\Python311\Lib\site-packages\pybind11\include\pybind11\numpy.h"


namespace py = pybind11;

// Average function using NumPy input
double cpp_average(py::array_t<double> input) {
    auto buf = input.unchecked<1>();  // Unsafe but fast access to 1D array
    double sum = 0.0;

    for (ssize_t i = 0; i < buf.shape(0); ++i)
        sum += buf(i);

    return sum / buf.shape(0);
}

// Python module definition
PYBIND11_MODULE(average, m) {
    m.def("cpp_average", &cpp_average, "Compute average of a NumPy array");
}
