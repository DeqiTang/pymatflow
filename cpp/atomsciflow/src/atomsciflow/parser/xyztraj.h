#ifndef atomsciflow_INCLUDE_atomsciflow_PARSER_XYZTRAJ_H_
#define atomsciflow_INCLUDE_atomsciflow_PARSER_XYZTRAJ_H_

#include <fstream>
#include <iostream>
#include <regex>
#include <cstdlib>
//#include <numbers>
#include <set>
#include <map>
#include <armadillo>
#include <iomanip>
#include <string>

#include "atomsciflow/base/crystal.h"

namespace atomsciflow {
//

int read_xyztraj_file(std::vector<atomsciflow::Crystal>& traj, std::string filepath, int readcell);


int write_xyztraj_file(std::vector<atomsciflow::Crystal>& traj, std::string filepath, int writecell);

} //end namespace atomsciflow


#endif // atomsciflow_INCLUDE_atomsciflow_PARSER_XYZTRAJ_H_


