
#include <iostream>
#include <pybind11/pybind11.h>

#include "atomsciflow/cp2k/cp2k.h"

namespace py = pybind11;


void add_class_cp2k(py::module& m) {
    py::class_<atomsciflow::Cp2k>(m, "Cp2k")
        .def(py::init())
        .def("to_string", &atomsciflow::Cp2k::to_string)
        //.def("set_subsys", &atomsciflow::Cp2k::set_subsys, py::return_value_policy::reference)
        ;
}

PYBIND11_MODULE(cp2k, m) {
    m.doc() = "cpptest module";
    m.attr("__version__") = "0.0.1";

    add_class_cp2k(m);
}
