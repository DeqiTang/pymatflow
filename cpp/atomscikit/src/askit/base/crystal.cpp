#include "askit/base/crystal.h"
#include "askit/base/atom.h"
#include "askit/base/element.h"
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



namespace askit {


int Crystal::read_xyz_file(std::string filepath) {
    
    int i = 0; // for iteration
    int j = 0; // for iteration
    
    int natom = 0; // number of atoms each image

    std::ifstream xyzfile;
    xyzfile.open(filepath);
    std::string line;
    std::vector<std::string> lines;

    while (std::getline(xyzfile, line)) {
        // ignore empty lines
        std::istringstream iss(line);
        std::vector<std::string> line_split((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
        if (line_split.size() != 0) {
            lines.push_back(line);
        }
    }
    xyzfile.close();
    
    std::regex whitespace("\\s+");
    
    std::vector<std::string> line_split(std::sregex_token_iterator(lines[0].begin(), lines[0].end(), whitespace, -1), 
        std::sregex_token_iterator());
    natom = int(std::atof(line_split[0].c_str()));
    
    if (natom != (lines.size() - 2)) {
        std::cout <<  "Warning!!!" << " number of atoms specified in the first line is not equal to total line number of file minus 2(empty lines are removed)" << std::endl;
        std::exit(1);
    }
    
    std::vector<std::string> cell_line_split(std::sregex_token_iterator(lines[1].begin(), lines[1].end(), whitespace, -1), 
        std::sregex_token_iterator());
        
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    a.push_back(std::atof(cell_line_split[1].c_str()));
    a.push_back(std::atof(cell_line_split[2].c_str()));
    a.push_back(std::atof(cell_line_split[3].c_str()));
    b.push_back(std::atof(cell_line_split[5].c_str()));
    b.push_back(std::atof(cell_line_split[6].c_str()));
    b.push_back(std::atof(cell_line_split[7].c_str()));
    c.push_back(std::atof(cell_line_split[9].c_str()));
    c.push_back(std::atof(cell_line_split[10].c_str()));
    c.push_back(std::atof(cell_line_split[11].c_str()));
    this->cell.push_back(a);
    this->cell.push_back(b);
    this->cell.push_back(c);        
    
            
    for (i = 0; i < natom; i++) {

        std::vector<std::string> line_split(std::sregex_token_iterator(lines[i+2].begin(), lines[i+2].end(), whitespace, -1), 
            std::sregex_token_iterator());
        
        Atom atom;
        
        //std::cout << line_split[0] << line_split[1] ;

        atom.set_name(line_split[0]);
        atom.set_x(std::atof(line_split[1].c_str()));
        atom.set_y(std::atof(line_split[2].c_str()));
        atom.set_z(std::atof(line_split[3].c_str()));
        this->atoms.push_back(atom);
    }    

    return 0;

}




int Crystal::write_xyz_file(std::string filepath) {
    std::ofstream xyzfile;
    xyzfile.open(filepath);

    xyzfile.setf(std::ios::fixed);
    xyzfile << this->atoms.size() << "\n";
    
    xyzfile << "cell: "
        << std::setprecision(9) << std::setw(15) << this->cell[0][0] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[0][1] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[0][2] << " | " 
        << std::setprecision(9) << std::setw(15) << this->cell[1][0] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[1][1] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[1][2] << " | " 
        << std::setprecision(9) << std::setw(15) << this->cell[2][0] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[2][1] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[2][2] << "\n";
    
    for (auto atom : this->atoms) {
        xyzfile <<// atom.name << " " << atom.x << " " << atom.y << " " << atom.z << "\n";
            atom.name << "\t"
            << std::setprecision(9) << std::setw(15) << atom.x << "\t"
            << std::setprecision(9) << std::setw(15) << atom.y << "\t"
            << std::setprecision(9) << std::setw(15) << atom.z << "\t" << "\n";
    }
    xyzfile.close();
    
    return 0;
}


int Crystal::read_xyz_str(std::string str) {
    
    int i = 0; // for iteration
    int j = 0; // for iteration
    
    int natom = 0; // number of atoms each image

    std::regex newline("\\n+");

    std::vector<std::string> lines(std::sregex_token_iterator(str.begin(), str.end(), newline, -1),
        std::sregex_token_iterator());

    
    std::regex whitespace("\\s+");
    
    std::vector<std::string> line_split(std::sregex_token_iterator(lines[0].begin(), lines[0].end(), whitespace, -1), 
        std::sregex_token_iterator());
    natom = int(std::atof(line_split[0].c_str()));
    
    if (natom != (lines.size() - 2)) {
        std::cout <<  "Warning!!!" << " number of atoms specified in the first line is not equal to total line number of file minus 2(empty lines are removed)" << std::endl;
        std::exit(1);
    }
    
    std::vector<std::string> cell_line_split(std::sregex_token_iterator(lines[1].begin(), lines[1].end(), whitespace, -1), 
        std::sregex_token_iterator());
        
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    a.push_back(std::atof(cell_line_split[1].c_str()));
    a.push_back(std::atof(cell_line_split[2].c_str()));
    a.push_back(std::atof(cell_line_split[3].c_str()));
    b.push_back(std::atof(cell_line_split[5].c_str()));
    b.push_back(std::atof(cell_line_split[6].c_str()));
    b.push_back(std::atof(cell_line_split[7].c_str()));
    c.push_back(std::atof(cell_line_split[9].c_str()));
    c.push_back(std::atof(cell_line_split[10].c_str()));
    c.push_back(std::atof(cell_line_split[11].c_str()));
    this->cell.push_back(a);
    this->cell.push_back(b);
    this->cell.push_back(c);        
    
            
    for (i = 0; i < natom; i++) {

        std::vector<std::string> line_split(std::sregex_token_iterator(lines[i+2].begin(), lines[i+2].end(), whitespace, -1), 
            std::sregex_token_iterator());
        
        Atom atom;
        
        //std::cout << line_split[0] << line_split[1] ;

        atom.set_name(line_split[0]);
        atom.set_x(std::atof(line_split[1].c_str()));
        atom.set_y(std::atof(line_split[2].c_str()));
        atom.set_z(std::atof(line_split[3].c_str()));
        this->atoms.push_back(atom);
    }    

    return 0;

}

int Crystal::write_xyz_str(std::string& str) { 
    str = this->write_xyz_str();
    return 0;
}

std::string Crystal::write_xyz_str() {
    std::stringstream xyz_str_stream;
    xyz_str_stream.setf(std::ios::fixed);

    xyz_str_stream << this->atoms.size() << "\n";
    
    xyz_str_stream << "cell: "
        << std::setprecision(9) << std::setw(15) << this->cell[0][0] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[0][1] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[0][2] << " | " 
        << std::setprecision(9) << std::setw(15) << this->cell[1][0] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[1][1] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[1][2] << " | " 
        << std::setprecision(9) << std::setw(15) << this->cell[2][0] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[2][1] << " " 
        << std::setprecision(9) << std::setw(15) << this->cell[2][2] << "\n";
    
    for (auto atom : this->atoms) {
        xyz_str_stream <<// atom.name << " " << atom.x << " " << atom.y << " " << atom.z << "\n";
            atom.name << "\t"
            << std::setprecision(9) << std::setw(15) << atom.x << "\t"
            << std::setprecision(9) << std::setw(15) << atom.y << "\t"
            << std::setprecision(9) << std::setw(15) << atom.z << "\t" << "\n";
    }
   
    return xyz_str_stream.str();
}


int Crystal::build_supercell(std::vector<int> n) {
    
    int i = 0;
    int j = 0;
    int k = 0;
    
    arma::mat new_cell(3, 3);
    
    new_cell.row(0) = n[0] * arma::conv_to<arma::rowvec>::from(this->cell[0]);
    new_cell.row(1) = n[1] * arma::conv_to<arma::rowvec>::from(this->cell[1]);
    new_cell.row(2) = n[2] * arma::conv_to<arma::rowvec>::from(this->cell[2]);
    
    //std::cout << "n1: " << n[0] << std::endl;
    //std::cout << "n2: " << n[1] << std::endl;
    //std::cout << "n3: " << n[2] << std::endl;
    
    // build supercell: replica in three vector one by one
    
    int natom_now = 0;
    double x, y, z;
    for (i = 0; i < 3; i++) {
        //
        natom_now = this->atoms.size();
        for (j = 0;j < (n[i] - 1); j++) {
            for (k = 0; k < natom_now; k++) {
                //std::cout << "testing" << std::endl;
                x = this->atoms[k].x + (j + 1) * this->cell[i][0]; // old cell at present
                y = this->atoms[k].y + (j + 1) * this->cell[i][1]; // old cell at present
                z = this->atoms[k].z + (j + 1) * this->cell[i][2]; // old cell at present
                //std::cout << "testing." << std::endl;
                askit::Atom atm;
                //std::cout << "testing.." << std::endl;
                atm.name = this->atoms[k].name;
                atm.x = x;
                atm.y = y;
                atm.z = z;
                //std::cout << "testing..." << std::endl;
                this->atoms.push_back(atm);
                //std::cout << "testing." << std::endl;
            }
        }
        
    }
    
    this->cell[0]  = arma::conv_to<std::vector<double>>::from(new_cell.row(0));
    this->cell[1]  = arma::conv_to<std::vector<double>>::from(new_cell.row(1));
    this->cell[2]  = arma::conv_to<std::vector<double>>::from(new_cell.row(2));
    
    return 0;
}

} //namespace askit
