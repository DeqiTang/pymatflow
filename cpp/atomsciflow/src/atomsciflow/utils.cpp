/*
providing tools for structure manipulation
*/

#include "atomsciflow/utils.h"

#include "atomsciflow/base/crystal.h"
#include "atomsciflow/base/atom.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <armadillo>

namespace atomsciflow {

std::string version() {
    return "0.0.0";
}

int move_along(atomsciflow::Crystal* structure, std::vector<int> atoms_to_move, std::vector<int> direc, double disp) {
    /*
    :param structure: an instance of atomsciflow::Crystal()
    :params atoms_to_move: a list of atoms to move, counting starts with 0
    :param direc: three int number to indicate the direction to move along
        namely the crystal orientation index
    :param disp: displacement of the atoms in unit of Angstrom
    */
    // do some checking
    if (*std::max_element(atoms_to_move.begin(), atoms_to_move.end()) > (structure->natom() - 1)) {
        std::cout << "===================================================================" << std::endl;
        std::cout << "                              WARNING \n";
        std::cout << "------------------------------------------------------------------------------\n";
        std::cout << "the atom you are trying to move is beyond the number of atoms in the structure\n";
        std::exit(1);
    }
    

    arma::vec direc_cartesian(3);
    arma::mat latcell(3, 3);
    std::vector<std::vector<double>> cell = structure->get_cell();
    latcell.row(0) = arma::conv_to<arma::rowvec>::from(cell[0]);
    latcell.row(1) = arma::conv_to<arma::rowvec>::from(cell[1]);
    latcell.row(2) = arma::conv_to<arma::rowvec>::from(cell[2]);
    direc_cartesian = latcell.row(0) * direc[0] + latcell.row(1) * direc[1] + latcell.row(2) * direc[2];
    
    // normalize
    double length = arma::norm(direc_cartesian);
    // length = np.sqrt(direc_cartesian[0]**2+direc_cartesian[1]**2+direc_cartesian[2]**2)
    direc_cartesian = direc_cartesian / length;

    double deltax = direc_cartesian[0] * disp;
    double deltay = direc_cartesian[1] * disp;
    double deltaz = direc_cartesian[2] * disp;
    
    for (auto i : atoms_to_move) {
        structure->atoms[i].x += deltax;
        structure->atoms[i].y += deltay;
        structure->atoms[i].z += deltaz;
    }
    return 0;
}

int remove_atoms(atomsciflow::Crystal* structure, std::vector<int> atoms_to_remove) {
    /*
    :param structure: an instance of atomsciflow::Crystal()
    :params atoms_to_remove: a list of atoms to remove, counting starts with 0
    */
    // do some checking
    if (*std::max_element(atoms_to_remove.begin(), atoms_to_remove.end()) > (structure->natom() - 1)) {
        std::cout << "===================================================================" << std::endl;
        std::cout << "                              WARNING \n";
        std::cout << "------------------------------------------------------------------------------\n";
        std::cout << "the atom you are trying to remove is beyond the number of atoms in the structure\n";
        std::exit(1);
    }
    
    structure->remove_atoms(atoms_to_remove);
    //end
    return 0;
}


int vacuum_layer(atomsciflow::Crystal* structure, int plane, double thickness) {
    /*
    :param structure: an instance of atomsciflow::Crystal()
    :param plane: the plane chosen to add vacuum layer, can be:
        1 -> ab plane
        2 -> ac plane
        3 -> bc plane
    :param thickness: thickness of the vacuum layer
    */
    // 
    //

    arma::mat latcell(3, 3);
    latcell.row(0) = arma::conv_to<arma::rowvec>::from(structure->cell[0]);
    latcell.row(1) = arma::conv_to<arma::rowvec>::from(structure->cell[1]);
    latcell.row(2) = arma::conv_to<arma::rowvec>::from(structure->cell[2]);    
    
    if (plane == 1) {
        std::cout << "the normal_of_ab" << std::endl;
        arma::vec normal_of_ab = arma::cross(arma::conv_to<arma::vec>::from(latcell.row(0)), arma::conv_to<arma::vec>::from(latcell.row(1)));
        std::cout << "got the normal_of_ab" << std::endl;
        // get the cosine of the angle between c and outer product of ab
        double cosangle = arma::dot(latcell.row(2), normal_of_ab) / arma::norm(latcell.row(2)) / arma::norm(normal_of_ab);
        double proj_c_on_ab_normal = arma::norm(latcell.row(2)) * std::abs(cosangle);
        std::vector<double> z_all;
        for (auto& atom : structure->atoms) {
            z_all.push_back(atom.z);
        }
        std::cout << "got all z for atoms" << std::endl;
        
        double factor = (*std::max_element(z_all.begin(), z_all.end()) - *std::min_element(z_all.begin(), z_all.end()) + thickness) / proj_c_on_ab_normal * 1.0;
        structure->cell[2] = arma::conv_to<std::vector<double>>::from(latcell.row(2) * factor);
    } else if (plane == 2) {
        arma::vec normal_of_ac = arma::cross(arma::conv_to<arma::vec>::from(latcell.row(0)), arma::conv_to<arma::vec>::from(latcell.row(2)));
        // get the cosine of the angle between b and outer product of ac
        double cosangle = arma::dot(latcell.row(1), normal_of_ac) / arma::norm(latcell.row(1)) / arma::norm(normal_of_ac);
        double proj_b_on_ac_normal = arma::norm(latcell.row(1)) * std::abs(cosangle);
        std::vector<double> y_all;
        for (auto& atom : structure->atoms) {
            y_all.push_back(atom.y);
        }
        double factor = (*std::max_element(y_all.begin(), y_all.end()) - *std::min_element(y_all.begin(), y_all.end()) + thickness) / proj_b_on_ac_normal * 1.0;
        structure->cell[1] = arma::conv_to<std::vector<double>>::from(latcell.row(1) * factor);
    } else if (plane == 3) {
        arma::vec normal_of_bc = arma::cross(arma::conv_to<arma::vec>::from(latcell.row(1)), arma::conv_to<arma::vec>::from(latcell.row(2)));
        // get the cosine of the angle between b and outer product of ac
        double cosangle = arma::dot(latcell.row(0), normal_of_bc) / arma::norm(latcell.row(0)) / arma::norm(normal_of_bc);
        double proj_a_on_bc_normal = arma::norm(latcell.row(0)) * std::abs(cosangle);
        std::vector<double> x_all;
        for (auto& atom : structure->atoms) {
            x_all.push_back(atom.x);
        }
        double factor = (*std::max_element(x_all.begin(), x_all.end()) - *std::min_element(x_all.begin(), x_all.end()) + thickness) / proj_a_on_bc_normal * 1.0;
        structure->cell[0] = arma::conv_to<std::vector<double>>::from(latcell.row(0) * factor);        
    }
    
    // end
    return 0;
}



int set_frac_min_to_zero(atomsciflow::Crystal* structure) {
    /*
    :return an object of atomsciflow::Crystal()
    Note:
        set the fractional coordinate minimum to zero, this is a way of standardize the cif file
    */
    
    int i = 0;

    // now calc the fractional coordinates
    std::vector<atomsciflow::Atom> atoms_frac{structure->atoms};
    
    arma::mat latcell(3, 3);
    latcell.row(0) = arma::conv_to<arma::rowvec>::from(structure->cell[0]);
    latcell.row(1) = arma::conv_to<arma::rowvec>::from(structure->cell[1]);
    latcell.row(2) = arma::conv_to<arma::rowvec>::from(structure->cell[2]);      
    
    
    
    arma::mat convmat = arma::inv(arma::trans(latcell));
    
    std::vector<double> cart_tmp(3);
    arma::vec frac_tmp(3);
    
    for (i = 0; i < structure->natom(); i++) {
        //
        cart_tmp.push_back(structure->atoms[i].x);
        cart_tmp.push_back(structure->atoms[i].y);
        cart_tmp.push_back(structure->atoms[i].z);
        //frac_tmp = arma::dot(convmat, arma::vec(cart_tmp));
        frac_tmp = convmat * arma::vec(cart_tmp);
        atoms_frac[i].x = frac_tmp[0];
        atoms_frac[i].y = frac_tmp[1];
        atoms_frac[i].z = frac_tmp[2];
    }

    
    // set the minimum of fractional coord to to 0
    std::vector<double> values;
    for (i = 0; i < structure->natom(); i++) {
        values.push_back(atoms_frac[i].x);
    }    
    double min_frac_x;
    min_frac_x = *std::min_element(values.begin(), values.end());
    
    values.clear();
    for (i = 0; i < structure->natom(); i++) {
        values.push_back(atoms_frac[i].y);
    }
    double min_frac_y = *std::min_element(values.begin(), values.end());
    
    values.clear();
    for (i = 0; i < structure->natom(); i++) {
        values.push_back(atoms_frac[i].z);
    }
    double min_frac_z = *std::min_element(values.begin(), values.end());
    
    
    for (i = 0; i < structure->natom(); i++) {
        atoms_frac[i].x -= min_frac_x;
        atoms_frac[i].y -= min_frac_y;
        atoms_frac[i].z -= min_frac_z;
    }
    
            
    // now convert coord of atom in atoms_frac_within_new_cell to cartesian
    
    
    // update convmat for frac to cart conversion
    convmat = arma::trans(latcell);
    
    arma::vec cart_final;
    for (i = 0; i < structure->natom(); i++) {
        //
        frac_tmp.at(0) = atoms_frac[i].x;
        frac_tmp.at(1) = atoms_frac[i].y;
        frac_tmp.at(2) = atoms_frac[i].z;
        //cart_final = arma::dot(convmat, frac_tmp);
        cart_final = convmat * frac_tmp;
        structure->atoms[i].set_x(cart_final[0]);
        structure->atoms[i].set_y(cart_final[1]);
        structure->atoms[i].set_z(cart_final[2]);
    }
    
    return 0;
    //
}
        
int set_frac_within_zero_and_one(atomsciflow::Crystal* structure) {
    /*
    :return an object of atomsciflow::Crystal()
    Note:
        set the fractional coordinate within the range of 0 and 1, this is a way of standardize the cif file
    */
    
    int i = 0; // for iteration
    int j = 0; // for iteration

    // now conv to the fractional coordinates
    arma::mat latcell(3, 3);
    latcell.row(0) = arma::conv_to<arma::rowvec>::from(structure->cell[0]);
    latcell.row(1) = arma::conv_to<arma::rowvec>::from(structure->cell[1]);
    latcell.row(2) = arma::conv_to<arma::rowvec>::from(structure->cell[2]);      
    
    arma::mat convmat = arma::inv(arma::trans(latcell));
    arma::vec tmp_arma_vec(3);
    //std::vector<double> tmp_cpp_vec(3);
    for (i = 0; i < structure->natom(); i ++) {
        tmp_arma_vec(0) = structure->atoms[i].x;
        tmp_arma_vec(1) = structure->atoms[i].y;
        tmp_arma_vec(2) = structure->atoms[i].z;
        //tmp_arma_vec = arma::dot(convmat, tmp_arma_vec);
        tmp_arma_vec = convmat * tmp_arma_vec;
        structure->atoms[i].x = tmp_arma_vec(0);
        structure->atoms[i].y = tmp_arma_vec(1);
        structure->atoms[i].z = tmp_arma_vec(2);
    }
    
    
    // set the fractional coordinates within 0 and 1
    for (i = 0; i < structure->atoms.size(); i++) {
        while (structure->atoms[i].x >= 1) {
            structure->atoms[i].x -= 1;
        }
        while (structure->atoms[i].x < 0) {
            structure->atoms[i].x += 1;
        }
        while (structure->atoms[i].y >= 1) {
            structure->atoms[i].y -= 1;
        }
        while (structure->atoms[i].y < 0) {
            structure->atoms[i].y += 1;
        }
        while (structure->atoms[i].z >= 1) {
            structure->atoms[i].z -= 1;
        }
        while (structure->atoms[i].z < 0) {
            structure->atoms[i].z += 1;
        }
    }


    // now convert coord of atom in atoms_frac_within_new_cell to cartesian
    // update convmat for convertion from frac to cart
    convmat = arma::trans(latcell);
    
    for (i = 0; i < structure->atoms.size(); i++) {
        tmp_arma_vec(0) = structure->atoms[i].x;
        tmp_arma_vec(1) = structure->atoms[i].y;
        tmp_arma_vec(2) = structure->atoms[i].z;
        //tmp_arma_vec = arma::dot(convmat, tmp_arma_vec);
        tmp_arma_vec = convmat * tmp_arma_vec;
        structure->atoms[i].x = tmp_arma_vec(0);
        structure->atoms[i].y = tmp_arma_vec(1);
        structure->atoms[i].z = tmp_arma_vec(2);
    }
    
    return 0;
}
    
    
int inverse_geo_center(atomsciflow::Crystal* structure) {
    /*
    calc the geometric center of the system and make an inversion against that center
    :param structure: an instance of atomsciflow::Crystal()
    */
    
    int i = 0; // for iteration
    
    // calc the geometric center
    double x = 0;
    double y = 0;
    double z = 0;
    
    for (auto atom : structure->atoms) {
        x += atom.x;
        y += atom.y;
        z += atom.z;
    }
    
    x /= structure->atoms.size();
    y /= structure->atoms.size();
    z /= structure->atoms.size();
    
    // now get the symmetry image against the geometric center
    for (i = 0; i < structure->natom(); i++) {
        structure->atoms[i].x = x * 2 - structure->atoms[i].x;
        structure->atoms[i].y = y * 2 - structure->atoms[i].y;
        structure->atoms[i].z = z * 2 - structure->atoms[i].z;
    }
    
    // end
    return 0;
}


int inverse_point(atomsciflow::Crystal* structure, std::vector<double> point) {
    /*
    calc the geometric center of the system and make an inversion against that center
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    :param point: the inverse center point, like [0.0, 0.0, 0.0]
    */
    
    
    // now get the symmetry image against the inverse center
    for (auto& atom : structure->atoms) {
        atom.x = point[0] * 2 - atom.x;
        atom.y = point[1] * 2 - atom.y;
        atom.z = point[2] * 2 - atom.z;
    }
    
    // end
    return 0;
}



int inverse_cell_center(atomsciflow::Crystal* structure) {
    /*
    make an inversion against the cell center
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    */
    
    
    int i = 0; // for iteration
    
    //first transfer to fractional coordinate and inverse against [0.5, 0.5, 0.5]
    arma::mat latcell(3, 3);
    latcell.row(0) = arma::conv_to<arma::rowvec>::from(structure->cell[0]);
    latcell.row(1) = arma::conv_to<arma::rowvec>::from(structure->cell[1]);
    latcell.row(2) = arma::conv_to<arma::rowvec>::from(structure->cell[2]);    
    
    arma::mat convmat = arma::inv(arma::trans(latcell));
    
    arma::vec arma_vec(3);
    
    for (i = 0; i < structure->natom(); i++) {
        arma_vec[0] = structure->atoms[i].x;
        arma_vec[1] = structure->atoms[i].y;
        arma_vec[2] = structure->atoms[i].z;
        //arma_vec = arma::dot(convmat, arma_vec);
        arma_vec = convmat * arma_vec;
        structure->atoms[i].x = arma_vec[0];
        structure->atoms[i].y = arma_vec[1];
        structure->atoms[i].z = arma_vec[2];
    }
    
    for (i = 0; i < structure->natom(); i++) {
        structure->atoms[i].x = 0.5 * 2 - structure->atoms[i].x;
        structure->atoms[i].y = 0.5 * 2 - structure->atoms[i].y;
        structure->atoms[i].z = 0.5 * 2 - structure->atoms[i].z;
    }
    
    // convert frac to cartesian again
    // update convmat for convertion from frac to cart
    convmat = arma::trans(latcell);
    
    for (i = 0; i < structure->natom(); i++) {
        arma_vec[0] = structure->atoms[i].x;
        arma_vec[1] = structure->atoms[i].y;
        arma_vec[2] = structure->atoms[i].z;
        //arma_vec = arma::dot(convmat, arma_vec);
        arma_vec = convmat * arma_vec;
        structure->atoms[i].x = arma_vec[0];
        structure->atoms[i].y = arma_vec[1];
        structure->atoms[i].z = arma_vec[2];
    }

    //
    return 0;
}





int rotate_along_axis(atomsciflow::Crystal* structure, std::vector<int> rotate_atoms, std::vector<int> axis) {
    /*
    rotate the specified atoms along the specified axis
    :param structure: an instance of pymatflow.structure.crystal.crystal()
    */
    return 0;
}
    
    
    
int enlarge_atoms(atomsciflow::Crystal* structure) {
    /*
    :return out:
        atoms: [
                ["C", 0.00000, 0.000000, 0.0000],
                ["O", 0.00000, 0.500000, 0.0000],
                ...
            ]
    Note: will enlarge the atoms in the unit cell along both a, b, c and -a, -b, -c direction.
        The goal is to make sure when the cell rotate in the 3D space, it will always be filled
        with atoms.
    */
    
    /*
    from pymatflow.base.atom import Atom
    #
    cell = copy.deepcopy(structure.cell)
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    
    n1 = np.ceil(np.max([a, b, c]) / a ) * 2 # maybe times 2 is not needed
    n2 = np.ceil(np.max([a, b, c]) / b ) * 2
    n3 = np.ceil(np.max([a, b, c]) / c ) * 2
    n = [int(n1), int(n2), int(n3)]
    print(n)
    
    atoms = copy.deepcopy(structure.atoms)
    # build supercell: replica in three vector one by one
    for i in range(3):
        natom_now = len(atoms)
        for j in range(n[i] - 1):
            for atom in atoms[:natom_now]:
                x = atom.x + float(j + 1) * structure.cell[i][0]
                y = atom.y + float(j + 1) * structure.cell[i][1]
                z = atom.z + float(j + 1) * structure.cell[i][2]
                atoms.append(Atom(atom.name, x, y, z))
        # replicate in the negative direction of structure.cell[i]
        for atom in atoms[:natom_now*n[i]]:
            x = atom.x - float(n[i]) * structure.cell[i][0]
            y = atom.y - float(n[i]) * structure.cell[i][1]
            z = atom.z - float(n[i]) * structure.cell[i][2]
            atoms.append(Atom(atom.name, x, y, z))
    return [[atom.name, atom.x, atom.y, atom.z] for atom in atoms]

    */
    return 0;
}



std::vector<atomsciflow::Atom> enlarge_atoms_new_cell(atomsciflow::Crystal* structure, std::vector<std::vector<double> > new_cell) {
    /*
    Note: will enlarge the atoms in the unit cell along both a, b, c and -a, -b, -c direction of the new_cell !!!
        but the cell is not redefined, the returned atoms is not used to form crystal, but to be 
        tailored by redefine_lattice function to get atoms for the redfined lattice.
        The goal is to make sure when the new cell rotate in the 3D space, it will always be filled
        with atoms.    
    */
    
    int i = 0;
    int j = 0;
    int k = 0;
    // std::cout << "testing..." << std::endl;
    //
    arma::mat latcell(3, 3);
    latcell.row(0) = arma::conv_to<arma::rowvec>::from(structure->cell[0]);
    latcell.row(1) = arma::conv_to<arma::rowvec>::from(structure->cell[1]);
    latcell.row(2) = arma::conv_to<arma::rowvec>::from(structure->cell[2]);  
    
    double a = arma::norm(latcell.row(0));
    double b = arma::norm(latcell.row(1));
    double c = arma::norm(latcell.row(2));
    
    arma::mat latcell_new(3, 3);
    latcell_new.row(0) = arma::conv_to<arma::rowvec>::from(new_cell[0]);
    latcell_new.row(1) = arma::conv_to<arma::rowvec>::from(new_cell[1]);
    latcell_new.row(2) = arma::conv_to<arma::rowvec>::from(new_cell[2]);      
    
    double new_a = arma::norm(latcell_new.row(0));
    double new_b = arma::norm(latcell_new.row(1));
    double new_c = arma::norm(latcell_new.row(2));


    std::cout << "new_a: " << new_a << "; new_b: " << new_b << "; new_c: " << new_c << std::endl;

    std::vector<double> cpp_vec(3);
    cpp_vec[0] = new_a;
    cpp_vec[1] = new_b;
    cpp_vec[2] = new_c;
    //std::cout << "testing..." << std::endl;
    int n1 = std::ceil(*std::max_element(cpp_vec.begin(), cpp_vec.end()) / a) * 2; // maye times 2 is not needed
    int n2 = std::ceil(*std::max_element(cpp_vec.begin(), cpp_vec.end()) / b) * 2;
    int n3 = std::ceil(*std::max_element(cpp_vec.begin(), cpp_vec.end()) / c) * 2;
    
    
    std::vector<int> n(3);
    n[0] = n1;
    n[1] = n2;
    n[2] = n3;
    
    
    std::vector<atomsciflow::Atom> atoms = structure->atoms;
    
    std::cout << "n1: " << n[0] << std::endl;
    std::cout << "n2: " << n[1] << std::endl;
    std::cout << "n3: " << n[2] << std::endl;
    
    //std::cout << atoms[0].name << std::endl;
    // build supercell: replica in three vector one by one
    std::cout << "atoms.size(): " << atoms.size() << std::endl;
    int natom_now = 0;
    double x, y, z;
    for (i = 0; i < 3; i++) {
        //
        natom_now = atoms.size();
        for (j = 0; j < (n[i] - 1); j++) {
            for (k = 0; k < natom_now; k++) {
                //std::cout << "testing" << std::endl;
                x = atoms[k].x + (j + 1) * structure->cell[i][0];
                y = atoms[k].y + (j + 1) * structure->cell[i][1];
                z = atoms[k].z + (j + 1) * structure->cell[i][2];
                //std::cout << "testing." << std::endl;
                atomsciflow::Atom atm;
                //std::cout << "testing.." << std::endl;
                atm.name = atoms[k].name;
                atm.x = x;
                atm.y = y;
                atm.z = z;
                //std::cout << "testing..." << std::endl;
                atoms.push_back(atm);
                //std::cout << "testing." << std::endl;
            }
        }

        //std::cout << "natom_now-> " << natom_now << std::endl;
        // std::cout << "n[" << i << "]-> " << n[i] << std::endl;
        // replicate in the negative direction of sturcture.cell[i]
        for (k = 0; k < (natom_now * n[i]); k++) {
            x = atoms[k].x - n[i] * structure->cell[i][0];
            y = atoms[k].y - n[i] * structure->cell[i][1];
            z = atoms[k].z - n[i] * structure->cell[i][2];
            atomsciflow::Atom atm;
            atm.name = atoms[k].name;
            atm.x = x;
            atm.y = y;
            atm.z = z;
            //std::cout << "testing.." << std::endl;
            atoms.push_back(atm);
        }
        
    }
    
    return atoms;
}
    

int redefine_lattice(atomsciflow::Crystal* structure, std::vector<int> a, std::vector<int> b, std::vector<int> c, double precision=1.0e-8) {
    /*
    :param a, b, c: new lattice vectors in terms of old.
        new_a = a[0] * old_a + a[1] * old_b + a[2] * old_c
        like a=[1, 0, 0], b=[0, 1, 0], c=[0, 0, 1] actually defines the
        same lattice as old.
    :param precision, a value that is less than 1 and infinitely close to 1
        used to judge whether one atom is in another periodic of the redefined cell

    Method:
        first make a large enough supercell, which guarantee that all the atoms in the new lattice are inside
        the supercell.
        then redfine the cell, and calc the fractional coord of all atoms with regarding the new cell
        finally remove those atoms who's fractional coord is not within range [0, 1), and we can convert fractional
        coords to cartesian.
    Note:
        relationship among convertion of coords. the most important point is that all coords actually have one common
        reference system, namely the General XYZ coordinate system. all the cell are defined with XYZ as reference,
        and the convmat build from the cell(with XYZ as reference) can be applied only to atoms also with XYZ as ference,
        finally we convert frac to cartesian using convmat also defined using cell with XYZ as reference, so we get 
        the cartesian with general XYZ as reference. In the last all the coord of atoms and cell have the general 
        XYZ  system as reference. So it works!
    */
    
    int i = 0; // for iteration
    int j = 0; // for iteration
    
    
    std::cout << "a: " << a[0] << " " << a[1] << " " << a[2] << std::endl;
    std::cout << "b: " << b[0] << " " << b[1] << " " << b[2] << std::endl;
    std::cout << "c: " << c[0] << " " << c[1] << " " << c[2] << std::endl;
    
    std::cout << "structure->cell[0]: " << structure->cell[0][0] << " " << structure->cell[0][1] << " " << structure->cell[0][2] << std::endl;
    std::cout << "structure->cell[1]: " << structure->cell[1][0] << " " << structure->cell[1][1] << " " << structure->cell[1][2] << std::endl;
    std::cout << "structure->cell[0]: " << structure->cell[2][0] << " " << structure->cell[2][1] << " " << structure->cell[2][2] << std::endl;
    
    arma::mat latcell_old(3, 3);
    latcell_old.row(0) = arma::conv_to<arma::rowvec>::from(structure->cell[0]);
    latcell_old.row(1) = arma::conv_to<arma::rowvec>::from(structure->cell[1]);
    latcell_old.row(2) = arma::conv_to<arma::rowvec>::from(structure->cell[2]);  
    
    
    std::cout << "latcell_old[0]: " << latcell_old(0, 0) << " " << latcell_old(0, 1) << " " << latcell_old(0, 2) << std::endl;
    std::cout << "latcell_old[1]: " << latcell_old(1, 0) << " " << latcell_old(1, 1) << " " << latcell_old(1, 2) << std::endl;
    std::cout << "latcell_old[2]: " << latcell_old(2, 0) << " " << latcell_old(2, 1) << " " << latcell_old(2, 2) << std::endl;    
    
    arma::mat latcell_new(3, 3);
    
    latcell_new.row(0) =  a[0] * latcell_old.row(0) + a[1] * latcell_old.row(1) + a[2] * latcell_old.row(2);
    latcell_new.row(1) =  b[0] * latcell_old.row(0) + b[1] * latcell_old.row(1) + b[2] * latcell_old.row(2);
    latcell_new.row(2) =  c[0] * latcell_old.row(0) + c[1] * latcell_old.row(1) + c[2] * latcell_old.row(2);
    
    std::cout << "latcell_new[0]: " << latcell_new(0, 0) << " " << latcell_new(0, 1) << " " << latcell_new(0, 2) << std::endl;
    std::cout << "latcell_new[1]: " << latcell_new(1, 0) << " " << latcell_new(1, 1) << " " << latcell_new(1, 2) << std::endl;
    std::cout << "latcell_new[2]: " << latcell_new(2, 0) << " " << latcell_new(2, 1) << " " << latcell_new(2, 2) << std::endl;        

    
    // enlarge the system
    std::vector<std::vector<double> > new_cell;
    new_cell.push_back(arma::conv_to<std::vector<double> >::from(latcell_new.row(0)));
    new_cell.push_back(arma::conv_to<std::vector<double> >::from(latcell_new.row(1)));
    new_cell.push_back(arma::conv_to<std::vector<double> >::from(latcell_new.row(2)));
    
    std::cout << "new_cell[0]: " << new_cell[0][0] << " " << new_cell[0][1] << " " << new_cell[0][2] << std::endl;
    std::cout << "new_cell[1]: " << new_cell[1][0] << " " << new_cell[1][1] << " " << new_cell[1][2] << std::endl;
    std::cout << "new_cell[2]: " << new_cell[2][0] << " " << new_cell[2][1] << " " << new_cell[2][2] << std::endl;
    
    //std::cout << "testing..." << std::endl;
    
    std::vector<atomsciflow::Atom> atoms_container = enlarge_atoms_new_cell(structure, new_cell);

    std::cout << "redefine_lattice: enlarge_atoms_new_cell finished!" << std::endl;
    
    // now calc the fractional coordinates of all atoms in atoms_container with new_cell as reference
    arma::mat convmat(3, 3);
    convmat = arma::inv(arma::trans(latcell_new));
    arma::vec arma_vec(3);
    
    for (i = 0; i < atoms_container.size(); i++) {
        arma_vec(0) = atoms_container[i].x;
        arma_vec(1) = atoms_container[i].y;
        arma_vec(2) = atoms_container[i].z;
        arma_vec = convmat * arma_vec;
        atoms_container[i].x = arma_vec(0);
        atoms_container[i].y = arma_vec(1);
        atoms_container[i].z = arma_vec(2);
    }
    
    
    std::cout << "redefine_lattice: atoms_container convert to frac coord finished!" << std::endl;
    std::cout << "atoms_container.size(): " << atoms_container.size() << std::endl;
    
    std::vector<atomsciflow::Atom> atoms_frac_within_new_cell;
    
    for (i = 0; i < atoms_container.size(); i++) {

        if ((0 <= atoms_container[i].x) && (atoms_container[i].x < (1 - precision))) {
            if ((0 <= atoms_container[i].y) && (atoms_container[i].y < (1 - precision))) {
                if ((0 <= atoms_container[i].z) && (atoms_container[i].z < (1 - precision))) {
                    atoms_frac_within_new_cell.push_back(atoms_container[i]);
                }
            }
        }
        
    }
    
    std::cout << "atoms_frac_within_new_cell->size(): " << atoms_frac_within_new_cell.size() << std::endl;
        
    // now convert coord of atom in atoms_frac_within_new_cell to cartesian
    
    convmat = arma::trans(latcell_new);
    structure->atoms.clear();
    
    for (i = 0; i < atoms_frac_within_new_cell.size(); i++) {
        arma_vec[0] = atoms_frac_within_new_cell[i].x;
        arma_vec[1] = atoms_frac_within_new_cell[i].y;
        arma_vec[2] = atoms_frac_within_new_cell[i].z;
        arma_vec = convmat * arma_vec;
        atomsciflow::Atom atm;
        atm.name = atoms_frac_within_new_cell[i].name;
        atm.x = arma_vec[0];
        atm.y = arma_vec[1];
        atm.z = arma_vec[2];
        structure->atoms.push_back(atm);
    }
    
    structure->cell = new_cell;
    
    //
    return 0;
}


    
int cleave_surface(atomsciflow::Crystal* structure, std::vector<int> direction, double thickness=10, double precision=1.0e-8) {
    /*
    :param structure: an instance of crystal()
    :param direction: direction of the surface plane, like [0, 0, 1], the reference of it is three lattice 
        vector a, b, c
    :param precision, a value that is less than 1 and infinitely close to 1
        used to judge whether one atom is in another periodic of the redefined cell used in cleave surface
        
    :return an object of crystal()
    
    Note: make use of redefine_lattice() to cleave surface.
        we try to find new_a and new_b which form the surface plane (the direction is normal to the plane)
    */
    
    arma::mat old_cell(3, 3);
    old_cell.row(0) = arma::conv_to<arma::rowvec>::from(structure->cell[0]);
    old_cell.row(1) = arma::conv_to<arma::rowvec>::from(structure->cell[1]);
    old_cell.row(2) = arma::conv_to<arma::rowvec>::from(structure->cell[2]);
    
    arma::vec direction_xyz_as_ref(3);
    direction_xyz_as_ref = arma::conv_to<arma::vec>::from(direction[0] * old_cell.row(0) + direction[1] * old_cell.row(1) + direction[2] * old_cell.row(2));
    std::cout << "direction_xyz_as_ref.size(): "  << direction_xyz_as_ref.size() << std::endl;        
    
    std::vector<int> a_from_old;
    std::vector<int> b_from_old;
    std::vector<int> c_from_old;

    std::vector<int> iter_ijk;
    
    int i = 0;
    int j = 0;
    int k = 0;
    
    iter_ijk.push_back(0);
    
    for (i = 0; i < 16; i++) {
        iter_ijk.push_back(i);
        iter_ijk.push_back(-i);
    }
    
    arma::vec new_vec(3);
    arma::vec new_a(3);
    
    std::cout << "begin search for new a and b" << std::endl;
    
    for (auto iter_i : iter_ijk) {
        i = iter_i;
        if ((! a_from_old.empty()) && (! b_from_old.empty())) {
            break;
        }
        
        for (auto iter_j : iter_ijk) {
            j = iter_j;
            if ((! a_from_old.empty()) && (! b_from_old.empty())) {
                break;
            }
            for (auto iter_k : iter_ijk) {
                k = iter_k;
                if ((! a_from_old.empty()) && (! b_from_old.empty())) {
                    break;
                }
                if ((i == 0) && (j == 0) && (k == 0 )) {
                    continue;
                }
                new_vec = arma::conv_to<arma::vec>::from(i * old_cell.row(0) + j * old_cell.row(1) + k * old_cell.row(2));
                //std::cout << "get new vec..." << std::endl;
                //std::cout << "new_vec.size(): "  << new_vec.size() << std::endl;
                //std::cout << "direction_xyz_as_ref.size(): "  << direction_xyz_as_ref.size() << std::endl;
                if (arma::dot(new_vec, direction_xyz_as_ref) == 0) {
                    if (a_from_old.empty()) {
                        a_from_old.push_back(i);
                        a_from_old.push_back(j);
                        a_from_old.push_back(k);
                        new_a = arma::conv_to<arma::vec>::from(a_from_old[0] * old_cell.row(0) + a_from_old[1] * old_cell.row(1) + a_from_old[2] * old_cell.row(2));
                        std::cout << "get new a..." << std::endl;
                        std::cout << "get a_from_old" << std::endl;
                        continue;
                    } else if (b_from_old.empty()) {
                        double cosangle = arma::dot(new_vec, new_a) / arma::norm(new_vec) / arma::norm(new_a);
                        double sinangle = std::sin(std::acos(cosangle));
                        std::cout << "sin: " << sinangle << std::endl;
                        std::cout << "cos: " << cosangle << std::endl;
                        if (std::abs(sinangle - 0) < 1.0e-3) {
                            continue;
                        } else if (1.0e-5 < cosangle && cosangle < 1) {
                            // ignore angle between [0, 90] by default
                            // like we ignore angle  = 60 degree for a and b here, we seek for 120 instead
                            continue;
                        } else {
                            std::cout << "get b_from_old" << std::endl;
                            b_from_old.push_back(i);
                            b_from_old.push_back(j);
                            b_from_old.push_back(k);
                            continue;
                        }
                    }
                }
            }
        }
    }
    
    
    std::cout << "begin search for new c" << std::endl;
    
    arma::vec new_c(3);
    std::vector<std::vector<int>> new_c_ijk;
    std::vector<double> new_c_sin;

    for (auto i_iter : iter_ijk) {
        i = i_iter;
        for (auto j_iter : iter_ijk) {
            j = j_iter;
            for (auto k_iter : iter_ijk) {
                k = k_iter;
                new_c = arma::conv_to<arma::vec>::from(i * old_cell.row(0) + j * old_cell.row(1) + k * old_cell.row(2));
                double cosangle = arma::dot(new_c, direction_xyz_as_ref) / arma::norm(new_c) / arma::norm(direction_xyz_as_ref);
                double sinangle = std::sin(std::acos(cosangle));
                std::vector<int> ijk;
                ijk.push_back(i);
                ijk.push_back(j);
                ijk.push_back(k);
                new_c_ijk.push_back(ijk);
                new_c_sin.push_back(sinangle);
            }
        }
    }
    
    std::vector<double> new_c_sin_remove_nan;
    for (i = 0; i < new_c_sin.size(); i++) {
        if (std::isnan(new_c_sin[i])) {
            continue;
        } else {
            new_c_sin_remove_nan.push_back(new_c_sin[i]);
        }
    }

    
    double min_sin = *std::min_element(new_c_sin_remove_nan.begin(), new_c_sin_remove_nan.end());
    
    std::cout << "min sin is: " << min_sin << std::endl;
    
    for (i = 0; i < new_c_sin.size(); i++) {
        if (new_c_sin[i] == min_sin) {
            c_from_old = new_c_ijk[i];
            break;
        }
    }


    std::cout << "a-> " << a_from_old[0] << " " << a_from_old[1] << " " << a_from_old[2] << std::endl;
    std::cout << "b-> " << b_from_old[0] << " " << b_from_old[1] << " " << b_from_old[2] << std::endl;
    std::cout << "c-> " << c_from_old[0] << " " << c_from_old[1] << " " << c_from_old[2] << std::endl;
    
    
    // redefine lattice
    std::cout << "redefine the lattice" << std::endl;
    redefine_lattice(structure, a_from_old, b_from_old, c_from_old, precision);
    
    std::cout << "add vacuum layer" << std::endl;
    vacuum_layer(structure, 1, thickness);
    
    return 0;
}




atomsciflow::Crystal merge_layers(atomsciflow::Crystal* structure1, atomsciflow::Crystal* structure2, int use_cell=0, double distance=3.4, double thickness=10) {
    /*
    :param structure1: an instance of crystal()
    :param structure2: an instance of crystal()    
    :param use_cell: use cell parameter of structure 1 or 2 or 0(average) to set the new a b cell parameter
        attention: c vector is not handled this way
        
    :param distance: the distance between layers
    :param thickness: the vaccum layer thickness of the combined system
    
    :return an object of crystal()
    Note:
        only merge layers with ab plane as the surface plane
    */
    
    int i = 0; // for iteration
    int j = 0; // for iteration
    int k = 0; // for iteration
    
    
    set_frac_within_zero_and_one(structure1);
    set_frac_within_zero_and_one(structure2);
    
    std::cout << "frac of structure 1 and 2 already set to be within zero and one" << std::endl;
    
    arma::mat old_cell_1(3, 3);
    old_cell_1.row(0) = arma::conv_to<arma::rowvec>::from(structure1->cell[0]);
    old_cell_1.row(1) = arma::conv_to<arma::rowvec>::from(structure1->cell[1]);
    old_cell_1.row(2) = arma::conv_to<arma::rowvec>::from(structure1->cell[2]);     
    
    arma::mat old_cell_2(3, 3);
    old_cell_2.row(0) = arma::conv_to<arma::rowvec>::from(structure2->cell[0]);
    old_cell_2.row(1) = arma::conv_to<arma::rowvec>::from(structure2->cell[1]);
    old_cell_2.row(2) = arma::conv_to<arma::rowvec>::from(structure2->cell[2]);     
    
    
    
    // first transfer to fractional coordinate
    std::cout << "first transfer to fractional coordinates" << std::endl;
    
    arma::mat latcell_old_1(3, 3);
    latcell_old_1.row(0) = arma::conv_to<arma::rowvec>::from(structure1->cell[0]);
    latcell_old_1.row(1) = arma::conv_to<arma::rowvec>::from(structure1->cell[1]);
    latcell_old_1.row(2) = arma::conv_to<arma::rowvec>::from(structure1->cell[2]);     
    
    arma::mat latcell_old_2(3, 3);
    latcell_old_2.row(0) = arma::conv_to<arma::rowvec>::from(structure2->cell[0]);
    latcell_old_2.row(1) = arma::conv_to<arma::rowvec>::from(structure2->cell[1]);
    latcell_old_2.row(2) = arma::conv_to<arma::rowvec>::from(structure2->cell[2]);     
    
    arma::mat convmat(3, 3);
    convmat = arma::inv(arma::trans(latcell_old_1));
    
    arma::vec arma_vec(3);
    
    std::vector<atomsciflow::Atom> frac_1;
    std::vector<atomsciflow::Atom> frac_2;
    
    for (i = 0; i < structure1->atoms.size(); i++) {
        arma_vec(0) = structure1->atoms[i].x;
        arma_vec(1) = structure1->atoms[i].y;
        arma_vec(2) = structure1->atoms[i].z;
        //arma_vec = arma::dot(convmat, arma_vec);
        arma_vec = convmat * arma_vec;
        atomsciflow::Atom atm;
        atm.name = structure1->atoms[i].name;
        atm.x = arma_vec(0);
        atm.y = arma_vec(1);
        atm.z = arma_vec(2);
        frac_1.push_back(atm);
    }
    
    convmat = arma::inv(arma::trans(latcell_old_2));
    
    for (i = 0; i < structure2->atoms.size(); i++) {
        arma_vec(0) = structure2->atoms[i].x;
        arma_vec(1) = structure2->atoms[i].y;
        arma_vec(2) = structure2->atoms[i].z;
        //arma_vec = arma::dot(convmat, arma_vec);
        arma_vec = convmat * arma_vec;
        atomsciflow::Atom atm;
        atm.name = structure2->atoms[i].name;
        atm.x = arma_vec(0);
        atm.y = arma_vec(1);
        atm.z = arma_vec(2);
        frac_2.push_back(atm);
    }

    //
    //
    arma::mat average_cell = (latcell_old_1 + latcell_old_2) / 2;
    
    
    atomsciflow::Crystal out;
    
    arma::mat latcell_frac_to_cart_1(3, 3);
    arma::mat latcell_frac_to_cart_2(3, 3);
    
    if (use_cell == 1) {
        latcell_frac_to_cart_1 = old_cell_1;
        latcell_frac_to_cart_2.row(0) = old_cell_1.row(0);
        latcell_frac_to_cart_2.row(1) = old_cell_1.row(1);
        latcell_frac_to_cart_2.row(2) = old_cell_2.row(2);
    } else if (use_cell == 2) {
        latcell_frac_to_cart_2 = old_cell_2;
        latcell_frac_to_cart_1.row(0) = old_cell_2.row(0);
        latcell_frac_to_cart_1.row(1) = old_cell_2.row(1);
        latcell_frac_to_cart_1.row(2) = old_cell_1.row(2);
    } else {
        //average_ab
        latcell_frac_to_cart_1.row(0) = (old_cell_1.row(0) + old_cell_2.row(0)) / 2;
        latcell_frac_to_cart_1.row(1) = (old_cell_1.row(1) + old_cell_2.row(1)) / 2;
        latcell_frac_to_cart_1.row(2) = old_cell_1.row(2);
        latcell_frac_to_cart_2.row(0) = (old_cell_1.row(0) + old_cell_2.row(0)) / 2;
        latcell_frac_to_cart_2.row(1) = (old_cell_1.row(1) + old_cell_2.row(1)) / 2;
        latcell_frac_to_cart_2.row(2) = old_cell_2.row(2);
    }
        

    // convert frac to cartesian again
    std::cout << "then transfer back to cartesian coordinates" << std::endl;
    arma::mat convmat_1 = arma::trans(latcell_frac_to_cart_1);
    arma::mat convmat_2 = arma::trans(latcell_frac_to_cart_2);
    
    std::vector<atomsciflow::Atom> cart_1;
    for (i = 0; i < frac_1.size(); i++) {
        arma_vec(0) = frac_1[i].x;
        arma_vec(1) = frac_1[i].y;
        arma_vec(2) = frac_1[i].z;
        //arma_vec = arma::dot(convmat_1, arma_vec);
        arma_vec = convmat_1 * arma_vec;
        atomsciflow::Atom atm;
        atm.name = frac_1[i].name;
        atm.x = arma_vec(0);
        atm.y = arma_vec(1);
        atm.z = arma_vec(2);
        cart_1.push_back(atm);
    }
    
    std::vector<atomsciflow::Atom> cart_2;
    for (i = 0; i < frac_2.size(); i++) {
        arma_vec(0) = frac_2[i].x;
        arma_vec(1) = frac_2[i].y;
        arma_vec(2) = frac_2[i].z;
        //arma_vec = arma::dot(convmat_2, arma_vec);
        arma_vec = convmat_2 * arma_vec;
        atomsciflow::Atom atm;
        atm.name = frac_2[i].name;
        atm.x = arma_vec(0);
        atm.y = arma_vec(1);
        atm.z = arma_vec(2);
        cart_2.push_back(atm);
    }
    
    // make distance gap between cart_1 and cart_2 is the value of distance
    std::cout << "make distance gap between cart_1 and cart_2 is the value of distance" << std::endl;
    std::vector<double> z_1;
    for (i = 0; i < cart_1.size(); i++) {
        z_1.push_back(cart_1[i].z);
    }
    
    std::vector<double> z_2;
    for (i = 0; i < cart_2.size(); i++) {
        z_2.push_back(cart_2[i].z);
    }
    
    double max_z_1 = *std::max_element(z_1.begin(), z_1.end());
    double min_z_2 = *std::min_element(z_2.begin(), z_2.end());
    
    std::cout << "max_z_1: " << max_z_1 << std::endl;
    std::cout << "min_z_1: " << min_z_2 << std::endl;
    
    for (i = 0; i < cart_2.size(); i++) {
        cart_2[i].z += distance - (min_z_2 - max_z_1);
    }

    
    //out.atoms = cart_1 + cart_2;
    std::cout << "merge atoms" << std::endl;
    out.atoms = cart_1;
    out.atoms.insert(out.atoms.end(), cart_2.begin(), cart_2.end());
    
    
    std::cout << "update cell information" << std::endl;
    
    double factor;
    
    if (use_cell == 1) {
        out.cell.clear();
        out.cell.push_back(arma::conv_to<std::vector<double>>::from(old_cell_1.row(0)));
        out.cell.push_back(arma::conv_to<std::vector<double>>::from(old_cell_1.row(1)));
        factor = (arma::norm(old_cell_1.row(2)) + arma::norm(old_cell_2.row(2))) / arma::norm(old_cell_1.row(2));
        out.cell.push_back(arma::conv_to<std::vector<double>>::from(old_cell_1.row(2) * factor));
    }  else if (use_cell == 2) {
        //
        out.cell.clear();
        out.cell.push_back(arma::conv_to<std::vector<double>>::from(old_cell_2.row(0)));
        out.cell.push_back(arma::conv_to<std::vector<double>>::from(old_cell_2.row(1)));
        factor = (arma::norm(old_cell_1.row(2)) + arma::norm(old_cell_2.row(2))) / arma::norm(old_cell_2.row(2));
        out.cell.push_back(arma::conv_to<std::vector<double>>::from(old_cell_2.row(2) * factor));        
    } else {
        // average cell
        out.cell.clear();
        out.cell.push_back(arma::conv_to<std::vector<double>>::from((old_cell_1.row(0) + old_cell_2.row(0)) / 2));
        out.cell.push_back(arma::conv_to<std::vector<double>>::from((old_cell_1.row(1) + old_cell_2.row(1)) / 2));
        factor = (arma::norm(old_cell_1.row(2)) + arma::norm(old_cell_2.row(2))) / arma::norm((old_cell_1.row(2) + old_cell_2.row(2)) / 2);
        out.cell.push_back(arma::conv_to<std::vector<double>>::from((old_cell_1.row(2) + old_cell_2.row(2)) / 2 * factor));
    }
    

    std::cout << "add vacuum layer: " << thickness << " (Angstrom)" << std::endl;
    vacuum_layer(&out, 1, thickness);
    std::cout << "end add vacuum layer" << std::endl;
    
    return out;
}


    
atomsciflow::Crystal build_nanotube_ab(atomsciflow::Crystal* structure, std::string axis="b") {
    /*
    :param structure: an instance of crystal()
    :param axis: a or b
    :return an object of crystal()
    Note:
        build nanotube along an axis parallel to axis a or b
        only apply to film structure with ab as plane and a must be vertical to b.
        ab plane must be periodical plane
    */
    
    int i = 0; // for iteration
    int j = 0; // for iteration
    int k = 0; // for iteration
    
    set_frac_within_zero_and_one(structure);
    
    atomsciflow::Crystal out;
    out.atoms = structure->atoms;
    out.cell = structure->cell;
    
    double a = arma::norm(arma::conv_to<arma::vec>::from(out.cell[0]));
    double b = arma::norm(arma::conv_to<arma::vec>::from(out.cell[1]));
    double c = arma::norm(arma::conv_to<arma::vec>::from(out.cell[2]));
    
    std::vector<double> all_x;
    std::vector<double> all_y;
    std::vector<double> all_z;

    for (i = 0; i < out.atoms.size(); i++) {
        all_x.push_back(out.atoms[i].x);
        all_y.push_back(out.atoms[i].y);
        all_z.push_back(out.atoms[i].z);
    }
    
    double radius;
    
    if (axis == "b") {
        radius = a / (2 * arma::datum::pi);
    } else if (axis == "a") {
        radius = b / (2 * arma::datum::pi);
    }
    
        
    double middle_z = *std::min_element(all_z.begin(), all_z.end()) + (*std::max_element(all_z.begin(), all_z.end()) - *std::min_element(all_z.begin(), all_z.end())) / 2;
    double center_x;
    double center_y;
    double center_z;
    if (axis == "b") {
        center_x = *std::min_element(all_x.begin(), all_x.end()) = a / 2;
        center_z = middle_z + radius;
    } else if (axis == "a") {
        center_y = *std::min_element(all_y.begin(), all_y.end()) + b / 2;
        center_z = middle_z + radius;
    }

    double x, y, z;

    for (i = 0; i < out.atoms.size(); i++) {
        x = out.atoms[i].x;
        y = out.atoms[i].y;
        z = out.atoms[i].z;
        if (axis == "b") {
            //
            double new_x, new_z;
            if (x >= center_x) {
                double arc_len = x - center_x;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_x = center_x + radius * std::sin(arc);
                new_z = center_z - radius * std::cos(arc);
                new_x += (middle_z - z) * std::cos(arc - arma::datum::pi / 2);
                new_z += (middle_z - z) * std::cos(arma::datum::pi - arc);
            }
            if (x < center_x) {
                double arc_len = center_x - x;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_x = center_x - radius * std::sin(arc);
                new_z = center_z - radius * std::cos(arc);
                new_x += (z - middle_z) * std::cos(arc - arma::datum::pi / 2);
                new_z -= (z - middle_z) * std::cos(arma::datum::pi - arc);
            }
            out.atoms[i].x = new_x;
            out.atoms[i].z = new_z;
        } else if (axis == "a") {
            //
            double new_y, new_z;
            if (y >= center_y) {
                double arc_len = y - center_y;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_y = center_y + radius * std::sin(arc);
                new_z = center_z - radius * std::cos(arc);
                new_y += (middle_z - z) * std::cos(arc - arma::datum::pi / 2);
                new_z += (middle_z - z) * std::cos(arma::datum::pi - arc);
            }
            if (y < center_y) {
                double arc_len = center_y - y;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_y = center_y - radius * std::sin(arc);
                new_z = center_z - radius * std::cos(arc);
                new_y += (z - middle_z) * std::cos(arc - arma::datum::pi / 2);
                new_z -= (z - middle_z) * std::cos(arma::datum::pi - arc);
            }
            out.atoms[i].y = new_y;
            out.atoms[i].z = new_z;
            
        }
    }
          
    double c_factor = 2 * radius / arma::norm(arma::conv_to<arma::vec>::from(out.cell[2])) * 2.5;
    out.cell[2] = arma::conv_to<std::vector<double>>::from(arma::conv_to<arma::vec>::from(out.cell[2]) * c_factor);
    
    return out;
}
    


atomsciflow::Crystal build_nanotube_ac(atomsciflow::Crystal* structure, std::string axis="a") {
    /*
    :param structure: an instance of crystal()
    :param axis: a or c
    :return an object of crystal()
    Note:
        build nanotube along an axis parallel to axis a or c
        only apply to film structure with ac as plane and a must be vertical to c.
        ac plane must be periodical plane        
    */
    
    int i = 0; // for iteration
    int j = 0; // for iteration
    
    atomsciflow::Crystal out;
    out.atoms = structure->atoms;
    out.cell = structure->cell;
    
    double a = arma::norm(arma::conv_to<arma::vec>::from(out.cell[0]));
    double b = arma::norm(arma::conv_to<arma::vec>::from(out.cell[1]));
    double c = arma::norm(arma::conv_to<arma::vec>::from(out.cell[2]));
    
    std::vector<double> all_x;
    std::vector<double> all_y;
    std::vector<double> all_z;
    
    for (i = 0; i < out.atoms.size(); i++) {
        all_x.push_back(out.atoms[i].x);
        all_y.push_back(out.atoms[i].y);
        all_z.push_back(out.atoms[i].z);
    }
    
    double radius;
    
    if (axis == "c") {
        radius = a / (2 * arma::datum::pi);
    } else if (axis == "a") {
        radius = c / (2 * arma::datum::pi);
    }
    
        
    double middle_y = *std::min_element(all_y.begin(), all_y.end()) + (*std::max_element(all_y.begin(), all_y.end()) - *std::min_element(all_y.begin(), all_y.end())) / 2;
    
    double center_x, center_y, center_z;
    
    if (axis == "c") {
        center_x = *std::min_element(all_x.begin(), all_x.end()) + a / 2;
        center_y = middle_y + radius;
    } else if (axis == "a") {
        center_z = *std::min_element(all_z.begin(), all_z.end()) + c / 2;
        center_y = middle_y + radius;
    }
    
    double x, y, z;
    
    for (i = 0; i < out.atoms.size(); i++) {
        x = out.atoms[i].x;    
        y = out.atoms[i].y;    
        z = out.atoms[i].z;
        if (axis == "c") {
            double new_x, new_y;
            if (x >= center_x) {
                double arc_len = x - center_x;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_x = center_x + radius * std::sin(arc);
                new_y = center_y - radius * std::cos(arc);
                new_x += (middle_y - y) * std::cos(arc - arma::datum::pi / 2);
                new_y += (middle_y - y) * std::cos(arma::datum::pi - arc);
            }
            if (x < center_x) {
                double arc_len = center_x - x;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_x = center_x - radius * std::sin(arc);
                new_y = center_y - radius * std::cos(arc);
                new_x += (y - middle_y) * std::cos(arc - arma::datum::pi / 2);
                new_y -= (y - middle_y) * std::cos(arma::datum::pi - arc);
            }
            out.atoms[i].x = new_x;
            out.atoms[i].y = new_y;
        } else if (axis == "a") {
            double new_y, new_z;
            if (z >= center_z) {
                double arc_len = z - center_z;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_z = center_z + radius * std::sin(arc);
                new_y = center_y - radius * std::cos(arc);
                new_z += (middle_y - y) * std::cos(arc - arma::datum::pi / 2);
                new_y += (middle_y - y) * std::cos(arma::datum::pi - arc);
            }
            if (z < center_z) {
                double arc_len = center_z - z;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_z = center_z - radius * std::sin(arc);
                new_y = center_y - radius * std::cos(arc);
                new_z += (y - middle_y) * std::cos(arc - arma::datum::pi / 2);
                new_y -= (y - middle_y) * std::cos(arma::datum::pi - arc);
            }
            out.atoms[i].y = new_y;
            out.atoms[i].z = new_z;
        }
    }
     
    double b_factor = 2 * radius / arma::norm(arma::conv_to<arma::vec>::from(out.cell[1])) * 2.5;
    out.cell[1] = arma::conv_to<std::vector<double>>::from(arma::conv_to<arma::vec>::from(out.cell[1]) * b_factor);
    
    return out;
}
    


atomsciflow::Crystal build_nanotube_bc(atomsciflow::Crystal* structure, std::string axis="c") {
    /*
    :param structure: an instance of crystal()
    :param axis: b or c
    :return an object of crystal()
    Note:
        build nanotube along an axis parallel to axis b or c
        only apply to film structure with bc as plane and b must be vertical to c.
        bc plane must be periodical plane        
    */
    
    
    int i = 0; // for iteration
    int j = 0; // for iteartion

    atomsciflow::Crystal out;
        
    out.atoms = structure->atoms;
    out.cell = structure->cell;
    
    double a = arma::norm(arma::conv_to<arma::vec>::from(out.cell[0]));
    double b = arma::norm(arma::conv_to<arma::vec>::from(out.cell[1]));
    double c = arma::norm(arma::conv_to<arma::vec>::from(out.cell[2]));
    
    std::vector<double> all_x;
    std::vector<double> all_y;
    std::vector<double> all_z;

    for (i = 0; i < out.atoms.size(); i++) {
        all_x.push_back(out.atoms[i].x);
        all_y.push_back(out.atoms[i].y);
        all_z.push_back(out.atoms[i].z);
    }
    
    double radius;
    if (axis == "b") {
        radius = c / (2 * arma::datum::pi);
    } else if (axis == "c") {
        radius = b / (2 * arma::datum::pi);
    }
    
        
    double middle_x = *std::min_element(all_x.begin(), all_x.end()) + (*std::max_element(all_x.begin(), all_x.end()) - *std::min_element(all_x.begin(), all_x.end())) / 2;
    
    double center_x, center_y, center_z;
    if (axis == "b") {
        center_z = *std::min_element(all_z.begin(), all_z.end()) + c / 2;
        center_x = middle_x + radius;
    } else if (axis == "c") {
        center_y = *std::min_element(all_y.begin(), all_y.end()) + b / 2;
        center_x = middle_x  + radius;
    }
    
    double x, y, z;
    
    for (i = 0; i < out.atoms.size(); i++) {
        x = out.atoms[i].x;
        y = out.atoms[i].y;
        z = out.atoms[i].z;
        if (axis == "b") {
            double new_x, new_z;
            if (z >= center_z) {
                double arc_len = z - center_z;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_z = center_z + radius * std::sin(arc);
                new_x = center_x - radius * std::cos(arc);
                new_z += (middle_x - x) * std::cos(arc - arma::datum::pi / 2);
                new_x += (middle_x - x) * std::cos(arma::datum::pi - arc);
            }
            if (z < center_z) {
                double arc_len = center_z - z;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_z = center_z - radius * std::sin(arc);
                new_x = center_x - radius * std::cos(arc);
                new_z += (x - middle_x) * std::cos(arc - arma::datum::pi / 2);
                new_x -= (x - middle_x) * std::cos(arma::datum::pi - arc);
            }
            out.atoms[i].x = new_x;
            out.atoms[i].z = new_z;
        } else if (axis == "c") {
            double new_x, new_y;
            if (y >= center_y) {
                double arc_len  = y - center_y;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_y = center_y + radius * std::sin(arc);
                new_x = center_x = radius * std::cos(arc);
                new_y += (middle_x - x) * std::cos(arc - arma::datum::pi / 2);
                new_x += (middle_x - x) * std::cos(arma::datum::pi - arc);
            }
            if (y < center_y) {
                double arc_len = center_y - y;
                // arc_len = radius * arc
                double arc = arc_len / radius;
                new_y = center_y - radius * std::sin(arc);
                new_x = center_x - radius * std::cos(arc);
                new_y += (x - middle_x) * std::cos(arc - arma::datum::pi / 2);
                new_x -= (x - middle_x) * std::cos(arma::datum::pi - arc);
            }
            out.atoms[i].y = new_y;
            out.atoms[i].x = new_x;
        }
    }
    
    
    double a_factor = 2 * radius / arma::norm(arma::conv_to<arma::vec>::from(out.cell[0])) * 2.5;
    out.cell[0] = arma::conv_to<std::vector<double>>::from(arma::conv_to<arma::vec>::from(out.cell[0]) * a_factor);
    
    return out;
}

} //namespace atomsciflow

