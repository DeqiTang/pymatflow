#ifndef ATOMSCIKIT_INCLUDE_ASKIT_PARSER_XYZ_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_PARSER_XYZ_H_

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
#include "askit/base/element.h"
#include "askit/base/atom.h"


namespace askit {
//

int read_xyz_file(askit::Crystal* crystal, std::string filepath);

int write_xyz_file(askit::Crystal* crystal, std::string filepath);

//int read_xyz_string(askit::Crystal* crystal, std::string str);

//int write_xyz_string(askit::Crystal* crystal, std::string& str);

//askit::Crystal read_xyz_string(std::string str);

//std::string write_xyz_string(askit::Crystal* crystal);

} // namespace askit

#endif // ATOMSCIKIT_INCLUDE_ASKIT_PASER_XYZ_H_
