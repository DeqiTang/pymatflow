#include "cpptest/test.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int a = 1, int b = 1) {
    return a + b;
}

PYBIND11_MODULE(pyaskit, m) {
    m.doc() = "cpptest function";
    m.def("add", &add, "a function adding to numbers", py::arg("i") = 1, py::arg("j") = 1);

    // export variables
    m.attr("the_answer") = 42;
    py::object world = py::cast("World");
    m.attr("what") = world;
}