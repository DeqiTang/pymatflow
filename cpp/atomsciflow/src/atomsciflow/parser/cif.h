#ifndef ATOMSCIFLOW_PARSER_CIF_H_
#define ATOMSCIFLOW_PARSER_CIF_H_

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

int read_cif_file(atomsciflow::Crystal* crystal, std::string filepath);

int write_cif_file(atomsciflow::Crystal* crystal, std::string filepath);


} // namespace atomsciflow

#endif // ATOMSCIFLOW_PASER_CIF_H_