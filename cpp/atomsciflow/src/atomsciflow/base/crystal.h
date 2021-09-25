/*
module for crystal structure manipulation
*/

#ifndef atomsciflow_INCLUDE_atomsciflow_BASE_CRYSTAL_H_
#define atomsciflow_INCLUDE_atomsciflow_BASE_CRYSTAL_H_

#include <vector>
#include <string>
#include <algorithm>
#include "atomsciflow/base/atom.h"

namespace atomsciflow {

class Crystal {

  public:

    Crystal() {};
    ~Crystal() {};
    
    int read_xyz_file(std::string filepath);

    int write_xyz_file(std::string filepath);

    int read_xyz_str(std::string str);

    int write_xyz_str(std::string& str);

    std::string write_xyz_str();

    std::vector<std::vector<double>> get_cell() {return this->cell;};
        /*
        ;params cell; {{a1, a2, a3}, {b1, b2, b3}, {c1, c2, c3}} in unit of Anstrom
        */

    int get_atoms(Atom* atoms) {return 0;};

    int get_cell_atoms(double** cell, Atom* atoms) {return 0;};


    int cartesian() {return 0;};

    int get_fractional() {return 0;};

    double volume() {return 0;};
    
    int build_supercell(std::vector<int> n);
    

    
    int to_base_xyz(std::string filepath) {return 0;};
    
    int remove_atom(int number) {return 0;};
   
    int remove_atoms(std::vector<int> numbers) {
      std::vector<Atom> origin = this->atoms;
      this->atoms.clear();
      
      for (int i = 0; i < origin.size(); i++) {
        //atoms_to_remove.push_back(this->atoms[i]);
        if (std::count(numbers.begin(), numbers.end(), i) == 0) {
          this->atoms.push_back(origin[i]);
        }
      }
      return 0;
    };

    int natom() { return this->atoms.size(); }
    
    std::vector<Atom> atoms;
    std::vector<std::vector<double>> cell;    
    
  private:


};

// cry.atoms[]

} //namespace atomsciflow

#endif // atomsciflow_INCLUDE_atomsciflow_BASE_CRYSTAL_H_
