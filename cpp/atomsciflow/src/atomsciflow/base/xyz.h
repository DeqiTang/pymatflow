/*
module for xyz structure
*/

#ifndef ATOMSCIFLOW_BASE_XYZ_H_
#define ATOMSCIFLOW_BASE_XYZ_H_

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include "atomsciflow/base/atom.h"

namespace atomsciflow {

class XYZ {

  public:

    XYZ() {};
    ~XYZ() {};
    
    int read_xyz_file(std::string filepath);

    int write_xyz_file(std::string filepath);

    void set_species_number();

    std::vector<std::vector<double>> get_cell() {return this->cell;};
        /*
        ;params cell; {{a1, a2, a3}, {b1, b2, b3}, {c1, c2, c3}} in unit of Anstrom
        */

    int get_atoms(Atom* atoms) {return 0;};

    int get_cell_atoms(double** cell, Atom* atoms) {return 0;};


    int cartesian() {return 0;};

    int get_fractional() {return 0;};

    double volume() {return 0;};

    int natom() { return this->atoms.size(); }
    
    std::string file;
    int nspecies;
    std::map<std::string, int> specie_labels;
    std::vector<Atom> atoms;
    std::vector<std::vector<double>> cell;    
    
  private:


};


} //namespace atomsciflow

#endif // ATOMSCIFLOW_BASE_XYZ_H_
