#pragma once
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>
#include <vector>

namespace py = pybind11;

class PythonObjectWrapper {
private:
    py::object instance;

    py::tuple to_tuple(const std::vector<py::object>& args) const {
        py::tuple tup(args.size());
        for (size_t i = 0; i < args.size(); ++i) {
            tup[i] = args[i];
        }
        return tup;
    }

public:
    PythonObjectWrapper(py::object py_class, const std::vector<py::object>& args = {}) {
        instance = py_class(*to_tuple(args));  // Unpack Python args correctly
    }

    py::object get_attr(const std::string& name) const {
        return instance.attr(name.c_str());
    }

    void set_attr(const std::string& name, py::object value) {
        instance.attr(name.c_str()) = value;
    }

    py::object call_method(const std::string& method_name, const std::vector<py::object>& args = {}) {
        return instance.attr(method_name.c_str())(*to_tuple(args));  // Proper call with unpacking
    }

    py::object raw() const {
        return instance;
    }
};
