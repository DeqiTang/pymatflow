
#include <iostream>
#include <pybind11/pybind11.h>

#include "cpptest/test.h"
#include "cpptest/test2.h"

namespace py = pybind11;

PYBIND11_MODULE(cpptest, m) {
    m.doc() = "cpptest module";
    m.def("add", &add, "a function adding to numbers", py::arg("i") = 1, py::arg("j") = 1);
    m.def("add2", &add, "a function adding to numbers", py::arg("i") = 1, py::arg("j") = 1);

    m.def("log1", 
        []() {
            return std::string("log1 function to test log string");
        }, 
        "a function to test log1 string"
    );
    m.def("log2", 
        []() {
            return log();
        }, 
        "a function to test log2 string"
    );


    // export variables
    m.attr("the_answer") = 42;
    py::object world = py::cast("World");
    m.attr("what") = world;
}
