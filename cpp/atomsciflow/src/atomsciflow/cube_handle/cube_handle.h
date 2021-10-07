#ifndef ATOMSCIFLOW_CUBE_HANDLE_H_
#define ATOMSCIFLOW_CUBE_HANDLE_H_

#include <iostream>
#include <vector>

#include "atomsciflow/parser/cube.h"

namespace atomsciflow {

int cube_diff_1d(std::vector<std::string> input_files, std::string output_file, std::vector<std::string> abscissa);

} // namespace atomsciflow

#endif // ATOMSCIFLOW_CUBE_HANDLE_H_