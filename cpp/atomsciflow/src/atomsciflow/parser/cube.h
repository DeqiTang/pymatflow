#ifndef ATOMSCIFLOW_PARSER_CUBE_H_
#define ATOMSCIFLOW_PARSER_CUBE_H_

#include <vector>
#include <string>
#include <algorithm>

#include <armadillo>

#include "atomsciflow/base/crystal.h"

namespace atomsciflow {

class CubeElectronDensity {
  
  public:
    
    void read_cube(std::string filepath);
    
    arma::cube data;

    atomsciflow::Crystal crystal;

    int ngridx, ngridy, ngridz;
};


} // namespace atomsciflow


#endif // ATOMSCIFLOW_PARSER_CUBE_H_