#ifndef atomsciflow_INCLUDE_atomsciflow_PARSER_XYZ_H_
#define atomsciflow_INCLUDE_atomsciflow_PARSER_XYZ_H_

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
#include "atomsciflow/base/element.h"
#include "atomsciflow/base/atom.h"


namespace atomsciflow {
//

int read_xyz_file(atomsciflow::Crystal* crystal, std::string filepath);

int write_xyz_file(atomsciflow::Crystal* crystal, std::string filepath);

//int read_xyz_string(atomsciflow::Crystal* crystal, std::string str);

//int write_xyz_string(atomsciflow::Crystal* crystal, std::string& str);

//atomsciflow::Crystal read_xyz_string(std::string str);

//std::string write_xyz_string(atomsciflow::Crystal* crystal);

} // namespace atomsciflow

#endif // atomsciflow_INCLUDE_atomsciflow_PASER_XYZ_H_
