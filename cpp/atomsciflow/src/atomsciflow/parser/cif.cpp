//

#include "atomsciflow/parser/cif.h"

#include "atomsciflow/base/atom.h"


namespace atomsciflow {

int read_cif_file(atomsciflow::Crystal* crystal, std::string filepath) {
    //
    
    std::ifstream ciffile;
    std::vector<std::string> lines;
    ciffile.open(filepath);
    
    std::string line;
    
    while (std::getline(ciffile, line)) {
        lines.push_back(line);
    }
    
    ciffile.close();
    
    // read the frac coords and a b c alpha beta gamma
    std::vector<Atom> atoms_frac;
    double a, b, c, alpha, beta, gamma;

    std::regex whitespace("\\s+");
    //std::regex whitespace(" ", std::regex::awk);
    
    int n_lines = lines.size();
    
    // for traverse
    int i = 0;
    
    int loop_of_atom_i = 0;
    int loop_atom_site_fract_x_i = 0;
    int loop_atom_site_fract_y_i = 0;
    int loop_atom_site_fract_z_i = 0;
    int col_name_i = 0;
    int col_x_i = 0;
    int col_y_i = 0;
    int col_z_i = 0;
    
    for (i = 0; i < lines.size(); i++) {

        //std::vector<std::string> line_split{std::sregex_token_iterator(lines[i].begin(), lines[i].end(), whitespace, -1), std::sregex_token_iterator{}};
        std::istringstream iss(lines[i]);
        std::vector<std::string> line_split((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
        
        if (line_split.size() == 0) {
            continue;
        }
        
        std::cout << "line->" << lines[i] << " line_split(size)->" << line_split.size() << "first->" << line_split[0] << std::endl;
        

        
        if (line_split[0] == "_cell_length_a") {
            a = std::atof(line_split[1].c_str());
            std::cout << "_cell_length_a: " << i << " value: " << a <<  std::endl;
            continue;
        } else if (line_split[0] == "_cell_length_b") {
            b = std::atof(line_split[1].c_str());
            continue;
        } else if (line_split[0] == "_cell_length_c") {
            c = std::atof(line_split[1].c_str());
            continue;
        } else if (line_split[0] == "_cell_angle_alpha") {
            alpha = std::atof(line_split[1].c_str());
            continue;
        } else if (line_split[0] == "_cell_angle_beta") {
            beta = std::atof(line_split[1].c_str());
            continue;
        } else if (line_split[0] == "_cell_angle_gamma") {
            gamma = std::atof(line_split[1].c_str());
            continue;
        } else if (lines[i].find("_atom_site_fract_x") != std::string::npos) {
            // not for sure about existing of \n
            std::cout << "_atom_site_fract_x: at line -> " << i << std::endl;
            loop_atom_site_fract_x_i = i;
            continue;
        } else if (lines[i].find("_atom_site_fract_y") != std::string::npos) {
            loop_atom_site_fract_y_i = i;
            continue;
        } else if (lines[i].find("_atom_site_fract_z") != std::string::npos) {
            loop_atom_site_fract_z_i = i;
            continue;
        }        
    }
    
    std::cout << "_atom_site_fract_x: " << loop_atom_site_fract_x_i << std::endl;
    std::cout << "_atom_site_fract_y: " << loop_atom_site_fract_y_i << std::endl;
    std::cout << "_atom_site_fract_z: " << loop_atom_site_fract_z_i << std::endl;
    
    // get loop_ line number for atom coords
    for (i = loop_atom_site_fract_x_i; i >=0; i--) {
        //std::vector<std::string> line_split(std::sregex_token_iterator(lines[i].begin(), lines[i].end(), whitespace, -1), std::sregex_token_iterator());
        if (lines[i].find("loop_") != std::string::npos) {
            loop_of_atom_i = i;
            break;
        }
    }
    
    // get ncol and begin of specific atom loop
    int ncol = 0;
    int loop_atom_block_truly_begin = 0;
    
    
    std::cout << "loop_of_atom_i: " << loop_of_atom_i << std::endl;
    
    for (i = (loop_of_atom_i + 1); i < n_lines; i++) {
        //std::vector<std::string> line_split{std::sregex_token_iterator(lines[i].begin(), lines[i].end(), whitespace, -1), std::sregex_token_iterator{}};
        std::istringstream iss(lines[i]);
        std::vector<std::string> line_split((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
        if (line_split.size() == 0) {
            continue;
        } else if (line_split.size() == 1) {
            ncol += 1;
            if (line_split[0] == "_atom_site_fract_x") {
                col_x_i = ncol - 1;
            } else if (line_split[0] == "_atom_site_fract_y") {
                col_y_i = ncol - 1;
            } else if (line_split[0] == "_atom_site_fract_z") {
                col_z_i = ncol - 1;
            } else if (line_split[0] == "_atom_site_type_symbol") {
                col_name_i = ncol - 1;
            }
        } else if (line_split.size() > 1) {
            loop_atom_block_truly_begin = i;
            break;
        }
    }
    
    
    std::cout << "col_x_i->" << col_x_i << "|" << "col_y_i->" << col_y_i << "|" << "col_z_i->" << col_z_i << "|" << "col_name_i->" << col_name_i << "|" << std::endl;
    
    std::string str_ = "_";
    
    std::cout << "loop_atom_block_truly_begin: " << loop_atom_block_truly_begin << std::endl;
    
    for (i = loop_atom_block_truly_begin; i < n_lines; i++) {
        //
        //std::vector<std::string> line_split{std::sregex_token_iterator(lines[i].begin(), lines[i].end(), whitespace, -1), std::sregex_token_iterator{}};
        std::istringstream iss(lines[i]);
        std::vector<std::string> line_split((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
        
        //std::cout << (line_split[0].find_first_of(str_) == 0) << std::endl;
        //std::cout << line_split[0] << std::endl;
        
        if (line_split.size() == 0) {
            continue;
        }
        
        if ((line_split.size() != ncol) || (line_split[0].find(str_) != std::string::npos)) {
            // loop alreay end
            break;
        }
        std::cout << line_split[0] << std::endl;
        Atom atom;
            
        //std::cout << line_split[0] << line_split[1] ;

        atom.set_name(line_split[col_name_i]);
        atom.set_x(std::atof(line_split[col_x_i].c_str()));
        atom.set_y(std::atof(line_split[col_y_i].c_str()));
        atom.set_z(std::atof(line_split[col_z_i].c_str()));
        atoms_frac.push_back(atom);        
    }
    

    std::cout << "cif file reading done" << std::endl;

    // convert alpha, beta, gamma to arc
    alpha = alpha  / 180 * arma::datum::pi;
    beta = beta / 180 * arma::datum::pi;
    gamma = gamma / 180 * arma::datum::pi;

    arma::mat latcell(3, 3);
    latcell(0, 0) = a;
    latcell(0, 1) = 0;
    latcell(0, 2) = 0;
    latcell(1, 0) = std::cos(gamma) * b;
    latcell(1, 1) = std::sin(gamma) * b;
    latcell(1, 2) = 0;
    double new_c1 = std::cos(beta);
    double new_c2 = (std::cos(alpha - std::cos(beta) * std::cos(gamma))) / std::sin(gamma);
    double new_c3_square = 1.0 - new_c1 * new_c1 - new_c2 * new_c2;
    double new_c3 = std::sqrt(new_c3_square);
    latcell(2, 0) = new_c1 * c;
    latcell(2, 1) = new_c2 * c;
    latcell(2, 2) = new_c3 * c;


    crystal->cell.clear();
    crystal->cell.push_back(arma::conv_to<std::vector<double>>::from(latcell.row(0)));
    crystal->cell.push_back(arma::conv_to<std::vector<double>>::from(latcell.row(1)));
    crystal->cell.push_back(arma::conv_to<std::vector<double>>::from(latcell.row(2)));
    
    // convert frac to cartesian
    std::cout << "converting frac to cartesian" << std::endl;
    arma::mat convmat = arma::trans(latcell);

    arma::vec xyz_cart(3);
    arma::vec xyz_frac(3);



    crystal->atoms.clear();
    for (auto& atom : atoms_frac) {
        Atom atm;
        xyz_frac[0] = atom.x;
        xyz_frac[1] = atom.y;
        xyz_frac[2] = atom.z;
        xyz_cart = convmat * xyz_frac;
        atm.set_name(atom.name);
        atm.set_xyz(
            xyz_cart[0],
            xyz_cart[1],
            xyz_cart[2]
        );
        crystal->atoms.push_back(atm);
    }
    
    
    std::cout << "cif parsing finished!" << std::endl;
    
    return 0;
}



int write_cif_file(atomsciflow::Crystal* crystal, std::string filepath) {
    std::ofstream ciffile;
    ciffile.open(filepath);
    // convert this->atoms to fractional coordinates
    // and standardize the cell with a transform matrix

    //arma::Mat<double> latcell =  arma::conv_to<arma::Mat<double>>::from(this->cell);
    arma::mat latcell(3, 3);
    latcell.row(0) = arma::conv_to<arma::rowvec>::from(crystal->cell[0]);
    latcell.row(1) = arma::conv_to<arma::rowvec>::from(crystal->cell[1]);
    latcell.row(2) = arma::conv_to<arma::rowvec>::from(crystal->cell[2]);

    arma::mat convmat = arma::inv(arma::trans(latcell));

    arma::vec xyz_cart(3);
    arma::vec xyz_frac(3);
    std::vector<Atom> atom_frac;

    for (auto& atom : crystal->atoms) {
        Atom atm;
        xyz_cart[0] = atom.x;
        xyz_cart[1] = atom.y;
        xyz_cart[2] = atom.z;
        xyz_frac = convmat * xyz_cart;
        atm.set_name(atom.name);
        atm.set_xyz(
            xyz_frac[0],
            xyz_frac[1],
            xyz_frac[2]
        );
        atom_frac.push_back(atm);
    }

    // frac already got, we need to put a along x to standardize the cell
    double a = std::sqrt(arma::accu(arma::square(latcell.row(0))));
    double b = std::sqrt(arma::accu(arma::square(latcell.row(1))));
    double c = std::sqrt(arma::accu(arma::square(latcell.row(2))));


    // angle in arc
    double alpha = std::acos(arma::dot(latcell.row(1), latcell.row(2)) / (b*c));
    double beta  = std::acos(arma::dot(latcell.row(0), latcell.row(2)) / (a*c));
    double gamma = std::acos(arma::dot(latcell.row(0), latcell.row(1)) / (a*b));

    // get new cell
    // new_a = a * [1, 0, 0]
    // new_b = b * [cos_gamma, sin_gamma, 0]
    // new_c1 = cos_beta
    // new_c2 = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    // new_c3_square = 1. - cx * cx - cy * cy must >= 0
    // new_c3 = sqrt(new_c3_square)
    // new_c = c * [new_c1, new_c2, new_c3]

    arma::mat newcell(3, 3);
    newcell(0, 0) = a;
    newcell(0, 1) = 0;
    newcell(0, 2) = 0;
    newcell(1, 0) = std::cos(gamma) * b;
    newcell(1, 1) = std::sin(gamma) * b;
    newcell(1, 2) = 0;
    double new_c1 = std::cos(beta);
    double new_c2 = (std::cos(alpha - std::cos(beta) * std::cos(gamma))) / std::sin(gamma);
    double new_c3_square = 1.0 - new_c1 * new_c1 - new_c2 * new_c2;
    double new_c3 = std::sqrt(new_c3_square);
    newcell(2, 0) = new_c1 * c;
    newcell(2, 1) = new_c2 * c;
    newcell(2, 2) = new_c3 * c;


    ciffile << "_cell_length_a" << "\t" << a << std::endl;
    ciffile << "_cell_length_b" << "\t" << b << std::endl;
    ciffile << "_cell_length_c" << "\t" << c << std::endl;
    ciffile << "_cell_angle_alpha" << "\t" << alpha * 180 / std::acos(-1)  << std::endl;
    ciffile << "_cell_angle_beta" << "\t" << beta * 180 / std::acos(-1) << std::endl;
    ciffile << "_cell_angle_gamma" << "\t" << gamma * 180 / std::acos(-1) << std::endl;
    ciffile << "\n" << std::endl;
    ciffile << "loop_" << std::endl;
    ciffile << "  _symmetry_equiv_pos_as_xyz" << std::endl;
    ciffile << "  \'x, y, z\'" << std::endl;
    ciffile << "\n" << std::endl;
    ciffile << "loop_" << std::endl;
    ciffile << "  _atom_site_label" << std::endl;
    ciffile << "  _atom_site_occupancy" << std::endl;
    ciffile << "  _atom_site_fract_x" << std::endl;
    ciffile << "  _atom_site_fract_y" << std::endl;
    ciffile << "  _atom_site_fract_z" << std::endl;
    ciffile << "  _atom_site_thermal_displace_type" << std::endl;
    ciffile << "  _atom_site_B_iso_or_equiv" << std::endl;
    ciffile << "  _atom_site_type_symbol" << std::endl;

    std::map<std::string, Element> ele_num_map = get_element_number_map();

    std::map<std::string, int> elements;
    for (auto atom : atom_frac) {
        elements[atom.name] = ele_num_map[atom.name].number;
    }

    std::set<std::pair<int, std::string>> elements_index;

    for (auto& atom : elements) {
        elements_index.emplace(atom.second, atom.first);
    }

    int i = 0;
    for (auto& elem : elements_index) {
        for (auto& atom : atom_frac) {
            if (atom.name == elem.second) {
                i += 1;
                ciffile << 
                    std::setw(5) << atom.name << i << "\t" <<
                    std::setw(6) << 1.000 << "\t" <<
                    std::setw(10) << atom.x << "\t" << 
                    std::setw(10) << atom.y << "\t" << 
                    std::setw(10) << atom.z << "\t" << 
                    "Biso" << "\t" << 1.000 << "\t" << 
                    std::setw(3) << atom.name << std::endl;
            }
        }
        i = 0;
    }

    //
    ciffile.close();
    return 0;
}



} // end namespace atomsciflow