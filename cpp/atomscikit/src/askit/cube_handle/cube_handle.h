#ifndef ATOMSCIKIT_INCLUDE_ASKIT_CUBE_HANDLE_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_CUBE_HANDLE_H_

#include <iostream>
#include <vector>

#include "askit/parser/cube.h"

namespace askit {

int cube_diff_1d(std::vector<std::string> input_files, std::string output_file, std::vector<std::string> abscissa);

} // namespace askit

#endif // ATOMSCIKIT_INCLUDE_ASKIT_CUBE_HANDLE_H_