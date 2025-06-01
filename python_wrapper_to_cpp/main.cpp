// main.cpp

#include <iostream>
#include "PythonObjectWrapper.h"

namespace py = pybind11;

int main() {
    py::scoped_interpreter guard{}; // Start the interpreter

    // Load the Python module and class
    py::module point3d_mod = py::module::import("point3d");
    py::object point3d_class = point3d_mod.attr("Point3D");

    // Create two Point3D Python objects via C++
    std::vector<py::object> args1 = {py::float_(1.0), py::float_(2.0), py::float_(3.0)};
    PythonObjectWrapper p1(point3d_class, args1);

    std::vector<py::object> args2 = {py::float_(4.0), py::float_(6.0), py::float_(8.0)};
    PythonObjectWrapper p2(point3d_class, args2);

    // Access attributes
    double x = p1.get_attr("x").cast<double>();
    double y = p1.get_attr("y").cast<double>();

    std::cout << "p1.x + p1.y = " << (x + y) << std::endl;

    // Call method
    py::object dist = p1.call_method("distance_to", {p2.raw()});
    std::cout << "Distance: " << dist.cast<double>() << std::endl;

    return 0;
}
