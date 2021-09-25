#include "atomsciflow/parser/cube.h"

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
#include <cmath>
#include <boost/algorithm/string.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "atomsciflow/base/crystal.h"
#include "atomsciflow/base/atom.h"
#include "atomsciflow/base/element.h"


namespace atomsciflow {

void CubeElectronDensity::read_cube(std::string filepath) {
    
    int i, j, k;

    std::ifstream cubefile;

    cubefile.open(filepath);

    std::string line;
    std::vector<std::string> cube;

    while (std::getline(cubefile, line)) {
        // ignore empty lines
        std::istringstream iss(line);
        std::vector<std::string> line_split((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
        if (line_split.size() != 0) {
            cube.push_back(line);
        }
    }
    cubefile.close();

    double bohr_to_angstrom = 0.529177249;
    // read structure info
    std::vector<std::string> vec_str;
    boost::split(vec_str, cube[2], boost::is_space(), boost::token_compress_on);
    int natom = int(std::atof(vec_str[1].c_str())); // it might be negative, if MO infor are included in cube file

    arma::mat latcell(3, 3);
    this->crystal.cell.clear();
    for ( i = 0; i < 3; i++) {
        boost::split(vec_str, cube[i+3], boost::is_space(), boost::token_compress_on);
        for (j = 0; j < 3; j++) {
            std::cout << vec_str[1] << " + " << vec_str[j+2] << " " << std::endl;
            latcell(i, j) = int(std::atof(vec_str[1].c_str())) * std::atof(vec_str[j+2].c_str()) * bohr_to_angstrom;
        }
        std::cout << "split string" << cube[i+3] << "\n";
        this->crystal.cell.push_back(arma::conv_to<std::vector<double>>::from(latcell.row(i)));
    }

    std::cout << this->crystal.cell[0][0] << " " << this->crystal.cell[0][1] << " " <<this->crystal.cell[0][2] << std::endl;
    std::cout << this->crystal.cell[1][0] << " " << this->crystal.cell[1][1] << " " <<this->crystal.cell[1][2] << std::endl;
    std::cout << this->crystal.cell[2][0] << " " << this->crystal.cell[2][1] << " " <<this->crystal.cell[2][2] << std::endl;

    std::string label;
    std::map<std::string, Element> ele_num_map = get_element_number_map();
    for (i = 0; i < natom; i++) {
        boost::split(vec_str, cube[i+6], boost::is_space(), boost::token_compress_on);
        int atomic_number = int(std::atof(vec_str[0].c_str()));
        for (auto e : ele_num_map) {
            if (ele_num_map[e.first].number == atomic_number) {
                label = e.first;
            }
        }
        Atom atm;
        atm.set_name(label);
        boost::split(vec_str, cube[i+6], boost::is_space(), boost::token_compress_on);
        atm.set_xyz(
            std::atof(vec_str[2].c_str()) * bohr_to_angstrom,
            std::atof(vec_str[3].c_str()) * bohr_to_angstrom,
            std::atof(vec_str[4].c_str()) * bohr_to_angstrom
        );
        this->crystal.atoms.push_back(atm);
    }

    double a = arma::norm(latcell.row(0), 2);
    double b = arma::norm(latcell.row(1), 2);
    double c = arma::norm(latcell.row(2), 2);
    
    // assume three cube file have the same ngridx ngridy and ngridz
    boost::split(vec_str, cube[3], boost::is_space(), boost::token_compress_on);
    int ngridx = int(std::atof(vec_str[1].c_str()));       
    boost::split(vec_str, cube[4], boost::is_space(), boost::token_compress_on);
    int ngridy = int(std::atof(vec_str[1].c_str()));       
    boost::split(vec_str, cube[5], boost::is_space(), boost::token_compress_on);
    int ngridz = int(std::atof(vec_str[1].c_str()));       
    
    this->ngridx = ngridx;
    this->ngridy = ngridy;
    this->ngridz = ngridz;
    
    std::cout << "ngrid: " << this->ngridx << " " << this->ngridy << " " << this->ngridz << std::endl;
    
    // read grid value
    std::string tmp;
    vec_str.clear();
    
    // std::cout << "natom: " << natom << std::endl;
    #pragma omp parallel for private(tmp)
    for (i = natom + 6; i < cube.size(); i++) {
        tmp = cube[i];
        // boost::algorithm::replace_all(tmp, "\n", ";");
        boost::algorithm::replace_all(tmp, "\n", " ");
        vec_str.push_back(tmp);
    }

    std::string tmp_str = boost::algorithm::join(vec_str, " ");


    arma::vec dat_vec(tmp_str);

    this->data.set_size(ngridx, ngridy, ngridz);
    // charge data in cube file is in shape (ngridx, ngridy, ngridz)
    // while charge in *CHG* file is in shape (ngzf, ngyf, ngxf)
    // they are different!        
    // grid data in cube is iterated in different compared to *CHG* of vasp
    // data = dat.reshape(ngridx, ngridy, ngridz);
    for (i = 0; i < ngridx; i++) {
        for (j = 0; j < ngridy; j++) {
            for (k = 0; k < ngridz; k++) {
                //this->data.at(i, j, k) = dat_vec.at(i*ngridx + j*ngridy + k*ngridz); // this is wrong
                this->data.at(i, j, k) = dat_vec.at(i*ngridy*ngridz + j*ngridx + k);
                // std::cout << this->data.at(i, j, k) << "\n";
            }
        }
    }
    
    std::cout << "data: last" << std::endl;
    std::cout << this->data.at(ngridx-1, ngridy-1, ngridz-1) << std::endl; 

}


} // namespace atomsciflow