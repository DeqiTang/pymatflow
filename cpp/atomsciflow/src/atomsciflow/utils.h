#ifndef ATOMSCIFLOW_UTILS_H_
#define ATOMSCIFLOW_UTILS_H_

#include "atomsciflow/base/crystal.h"
#include "atomsciflow/base/atom.h"


int move_along(atomsciflow::Crystal* structure, std::vector<int> atoms_to_move, std::vector<int> direc, double disp);

int remove_atoms(atomsciflow::Crystal* structure, std::vector<int> atoms_to_remove);

int vacuum_layer(atomsciflow::Crystal* structure, int plane, double thickness);

int set_frac_min_to_zero(atomsciflow::Crystal* structure);

int set_frac_within_zero_and_one(atomsciflow::Crystal* structure);

int inverse_geo_center(atomsciflow::Crystal* structure);

int inverse_point(atomsciflow::Crystal* structure, std::vector<double> point);

int inverse_cell_center(atomsciflow::Crystal* structure);

int rotate_along_axis(atomsciflow::Crystal* structure, std::vector<int> rotate_atoms, std::vector<int> axis);

int enlarge_atoms(atomsciflow::Crystal* structure);

std::vector<atomsciflow::Atom> enlarge_atoms_new_cell(atomsciflow::Crystal* structure, std::vector<std::vector<double> > new_cell);

int redefine_lattice(atomsciflow::Crystal* structure, std::vector<int> a, std::vector<int> b, std::vector<int> c, double precision);

int cleave_surface(atomsciflow::Crystal* structure, std::vector<int> direction, double thickness, double precision);

atomsciflow::Crystal merge_layers(atomsciflow::Crystal* structure1, atomsciflow::Crystal* structure2, int use_cell, double distance, double thickness);

atomsciflow::Crystal build_nanotube_ab(atomsciflow::Crystal* structure, std::string axis);

atomsciflow::Crystal build_nanotube_ac(atomsciflow::Crystal* structure, std::string axis);

atomsciflow::Crystal build_nanotube_bc(atomsciflow::Crystal* structure, std::string axis);


#endif // ATOMSCIFLOW_UTILS_H_