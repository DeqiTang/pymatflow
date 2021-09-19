#ifndef ATOMSCIKIT_INCLUDE_ASKIT_PARSER_CIF_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_PARSER_CIF_H_

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

int read_cif_file(askit::Crystal* crystal, std::string filepath);

int write_cif_file(askit::Crystal* crystal, std::string filepath);


} // namespace askit

#endif // ATOMSCIKIT_INCLUDE_ASKIT_PASER_CIF_H_