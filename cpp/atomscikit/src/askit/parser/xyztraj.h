#ifndef ATOMSCIKIT_INCLUDE_ASKIT_PARSER_XYZTRAJ_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_PARSER_XYZTRAJ_H_

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

#include "askit/base/crystal.h"

namespace askit {
//

int read_xyztraj_file(std::vector<askit::Crystal>& traj, std::string filepath, int readcell);


int write_xyztraj_file(std::vector<askit::Crystal>& traj, std::string filepath, int writecell);

} //end namespace askit


#endif // ATOMSCIKIT_INCLUDE_ASKIT_PARSER_XYZTRAJ_H_


