#ifndef ATOMSCIKIT_INCLUDE_ASKIT_UTILS_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_UTILS_H_

#include "askit/base/crystal.h"
#include "askit/base/atom.h"


int move_along(askit::Crystal* structure, std::vector<int> atoms_to_move, std::vector<int> direc, double disp);

int remove_atoms(askit::Crystal* structure, std::vector<int> atoms_to_remove);

int vacuum_layer(askit::Crystal* structure, int plane, double thickness);

int set_frac_min_to_zero(askit::Crystal* structure);

int set_frac_within_zero_and_one(askit::Crystal* structure);

int inverse_geo_center(askit::Crystal* structure);

int inverse_point(askit::Crystal* structure, std::vector<double> point);

int inverse_cell_center(askit::Crystal* structure);

int rotate_along_axis(askit::Crystal* structure, std::vector<int> rotate_atoms, std::vector<int> axis);

int enlarge_atoms(askit::Crystal* structure);

std::vector<askit::Atom> enlarge_atoms_new_cell(askit::Crystal* structure, std::vector<std::vector<double> > new_cell);

int redefine_lattice(askit::Crystal* structure, std::vector<int> a, std::vector<int> b, std::vector<int> c, double precision);

int cleave_surface(askit::Crystal* structure, std::vector<int> direction, double thickness, double precision);

askit::Crystal merge_layers(askit::Crystal* structure1, askit::Crystal* structure2, int use_cell, double distance, double thickness);

askit::Crystal build_nanotube_ab(askit::Crystal* structure, std::string axis);

askit::Crystal build_nanotube_ac(askit::Crystal* structure, std::string axis);

askit::Crystal build_nanotube_bc(askit::Crystal* structure, std::string axis);


#endif // ATOMSCIKIT_INCLUDE_ASKIT_UTILS_H_