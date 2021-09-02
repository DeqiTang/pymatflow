#include "askit/parser/xyztraj.h"

#include <iostream>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "askit/base/crystal.h"

namespace askit {

int read_xyztraj_file(std::vector<askit::Crystal>& traj, std::string filepath, int readcell) {
    /*
    cell must be specified in the xyz trajectory file
    :param readcell: 0 -> do not read cell; 1 -> read cell
    */
    
    traj.clear();
    
    int i = 0; // for iteration
    int j = 0; // for iteration
    
    int natom = 0; // number of atoms each image
    int nimage = 0; // number of images in the file
    
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


    //std::regex whitespace("\\s+");
    
    std::istringstream iss(lines[0]);
    std::vector<std::string> line_split((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
    natom = int(std::atof(line_split[0].c_str()));

    std::cout << "natom: " << natom << std::endl;
    std::cout << "number of lines(remove empty line): " << lines.size() << std::endl;
    nimage = int(lines.size() / (natom +2));

    std::cout << "num of images in xyz: " << nimage << std::endl;
    
    #pragma omp parallel for ordered
    for (i = 0; i < nimage; i++) {
        Crystal crystal;    
        if (readcell == 1) {
            //std::vector<std::string> line_split(std::sregex_token_iterator(lines[i*(natom+2)+1].begin(), lines[i*(natom+2)+1].end(), whitespace, -1), 
            //    std::sregex_token_iterator());    

            std::istringstream iss(lines[i*(natom+2)+1]);
            std::vector<std::string> line_split((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());


            std::vector<double> a;
            std::vector<double> b;
            std::vector<double> c;
            a.push_back(std::atof(line_split[1].c_str()));
            a.push_back(std::atof(line_split[2].c_str()));
            a.push_back(std::atof(line_split[3].c_str()));
            b.push_back(std::atof(line_split[5].c_str()));
            b.push_back(std::atof(line_split[6].c_str()));
            b.push_back(std::atof(line_split[7].c_str()));
            c.push_back(std::atof(line_split[9].c_str()));
            c.push_back(std::atof(line_split[10].c_str()));
            c.push_back(std::atof(line_split[11].c_str()));
            crystal.cell.push_back(a);
            crystal.cell.push_back(b);
            crystal.cell.push_back(c);
        }
        //
        for (j = 0; j < natom; j++) {
            Atom atom;
            
            //std::vector<std::string> line_split(std::sregex_token_iterator(lines[i*(natom+2)+2+j].begin(), lines[i*(natom+2)+2+j].end(), whitespace, -1), 
            //    std::sregex_token_iterator());    
            
            std::istringstream iss(lines[i*(natom+2)+2+j]);
            std::vector<std::string> line_split((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
            
            atom.set_name(line_split[0]);
            atom.set_x(std::atof(line_split[1].c_str()));
            atom.set_y(std::atof(line_split[2].c_str()));
            atom.set_z(std::atof(line_split[3].c_str()));
            crystal.atoms.push_back(atom);            
        }
        #pragma omp ordered
        traj.push_back(crystal);
    }


    return 0;
}



int write_xyztraj_file(std::vector<askit::Crystal>& traj, std::string filepath, int writecell) {
    //
    std::ofstream xyzfile;
    xyzfile.open(filepath);

    xyzfile.setf(std::ios::fixed);
    
    for (int i = 0; i < traj.size(); i++) {
        xyzfile << traj[i].atoms.size() << "\n";
    
        if (writecell == 1) {
            xyzfile << "cell: "
                << std::setprecision(9) << std::setw(15) << traj[i].cell[0][0] << " " 
                << std::setprecision(9) << std::setw(15) << traj[i].cell[0][1] << " " 
                << std::setprecision(9) << std::setw(15) << traj[i].cell[0][2] << " | " 
                << std::setprecision(9) << std::setw(15) << traj[i].cell[1][0] << " " 
                << std::setprecision(9) << std::setw(15) << traj[i].cell[1][1] << " " 
                << std::setprecision(9) << std::setw(15) << traj[i].cell[1][2] << " | " 
                << std::setprecision(9) << std::setw(15) << traj[i].cell[2][0] << " " 
                << std::setprecision(9) << std::setw(15) << traj[i].cell[2][1] << " " 
                << std::setprecision(9) << std::setw(15) << traj[i].cell[2][2] << "\n";
        } else {
            xyzfile << i << "-th image" << "\n";
        }

        for (auto& atom : traj[i].atoms) {
            xyzfile <<// atom.name << " " << atom.x << " " << atom.y << " " << atom.z << "\n";
                atom.name << "\t"
                << std::setprecision(9) << std::setw(15) << atom.x << "\t"
                << std::setprecision(9) << std::setw(15) << atom.y << "\t"
                << std::setprecision(9) << std::setw(15) << atom.z << "\t" << "\n";
        }

    }
    
    xyzfile.close();

    return 0;    
}


} // end namespace askit