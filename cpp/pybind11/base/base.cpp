
#include <iostream>
#include <pybind11/pybind11.h>

#include "cpptest/test.h"
#include "cpptest/test2.h"
#include "atomsciflow/base/crystal.h"

namespace py = pybind11;


void add_class_crystal(py::module& m) {
    py::class_<atomsciflow::Crystal>(m, "Crystal")
        .def(py::init())
        .def("read_xyz_str", &atomsciflow::Crystal::read_xyz_str, py::return_value_policy::reference)
        ;
}

PYBIND11_MODULE(base, m) {
    m.doc() = "base module";
    m.attr("__version__") = "0.0.1";

    add_class_crystal(m);
}
