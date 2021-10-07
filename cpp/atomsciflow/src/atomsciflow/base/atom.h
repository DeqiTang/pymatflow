/*
atoms are a collection of atom
*/

#ifndef ATOMSCIFLOW_BASE_ATOM_H_
#define ATOMSCIFLOW_BASE_ATOM_H_

//#include "element.h"
#include <vector>
#include <string>

namespace atomsciflow {

class Atom {
  public:
    
    explicit Atom() {};

    ~Atom() {};

    std::string get_name(std::string name) { return this->name; }

    double get_x() { return this->x; }

    double get_y() { return this->y; }

    double get_z() { return this->z; }

    void set_name(std::string name) { this->name = name; }

    void set_x(double x) { 
        //this->cartesian[0] = x; 
        this->x = x;
    }

    void set_y(double y) {
        this->y = y;
    }

    void set_z(double z) {
        this->z = z;
    }

    void set_xyz(double x, double y, double z) {
        //this->cartesian[0] = x;
        //this->cartesian[1] = y;
        //this->cartesian[2] = z;
        this->x = x;
        this->y = y;
        this->z = z;
    }

    double x, y, z;
    std::string name;
    
  private:
    //double cartesian[3];

};


} //namespace atomsciflow

#endif // ATOMSCIFLOW_BASE_ATOM_H_