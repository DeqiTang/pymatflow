#include <iostream>
#include <pybind11/pybind11.h>
// for automatic recognization an converssion of vector map arguments
#include <pybind11/stl.h>

#include "askit/cube_handle/cube_handle.h"

namespace py = pybind11;

PYBIND11_MODULE(cube_handle, m) {
    m.doc() = "cube handling module";
    m.def("cube_diff_1d", &askit::cube_diff_1d, py::return_value_policy::reference);
}