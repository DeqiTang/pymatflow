#ifndef ATOMSCIKIT_INCLUDE_ASKIT_PARSER_CUBE_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_PARSER_CUBE_H_

#include <vector>
#include <string>
#include <algorithm>

#include <armadillo>

#include "askit/base/crystal.h"

namespace askit {

class CubeElectronDensity {
  
  public:
    
    void read_cube(std::string filepath);
    
    arma::cube data;

    askit::Crystal crystal;

    int ngridx, ngridy, ngridz;
};


} // namespace askit


#endif // ATOMSCIKIT_INCLUDE_ASKIT_PARSER_CUBE_H_