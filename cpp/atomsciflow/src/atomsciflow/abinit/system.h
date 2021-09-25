/************************************************************************
    > File Name: system.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sun 24 Jan 2021 07:59:44 PM CST
************************************************************************/

#ifndef atomsciflow_INCLUDE_ASKIT_ABINIT_SYSTEM_H_
#define atomsciflow_INCLUDE_ASKIT_ABINIT_SYSTEM_H_

#include <string>
#include "atomsciflow/base/xyz.h"

namespace atomsciflow {

class AbinitSystem {
    /*
     */
    pulibc:

    AbinitSystem() {
        this->xyz = new XYZ();
        this->status = true;
        this->coordtype = "cartesian"; // can be cartesian or reduced
    };
    //
    //
    std::string to_string(int n = 0) {
        /*
        *return input_str is the string of all the set params
        *
        */
        std::string intput_str = "";
        input_str += "# ==================================\n";
        input_str += "# system related setting\n";
        input_str += "# ==================================\n";
        input_str += "\n";
        // scaling with 1 means no actually scaling of rprim by acell
        if (n > 0) {
            input_str += "acell" + std::to_string(n) + " 1 1 1 angstrom\n\n";
            input_str += "rprim" + std::to_string(n) + "\n";
            for (int i = 0; i < 3; i++) {
                input_str += std::to_string(cell[i][0]) + " " + std::to_string(cell[i][1]) + " " + std::to_string(cell[i][2]) +"\n";
            }
            input_str += "\n";
            input_str += "ntypat" + std::to_string(n) + " " + std::to_string(this->xyz.nspecies) + "\n\n";
            input_str += "natom" + std::to_string(n) + " " + std::to_string(this->xyz.natom()) + "\n\n";

            input_str += "typat" + std::to_string(n) + "\n";
            // abinit 不允许输入文件列数超过264, 因此如果原子数太多
            // 这里的typat要分多行列出
            // 利用余数设置如果一行超过30个原子就换行
            int i = 0;
            for (auto atom : this->xyz.atoms) {
                input_str += std::to_string(this->xyz.specie_labels[atom.name]) + " ";
                if (i % 30 == 29) {
                    input_str += "\n";
                }
                i += 1;
            }
            input_str += "\n\n";

            input_str += "znucl" + std::to_string(n) + "\n";

            for (auto element :this->xyz.specie_labels) {
                input_str += std::to_string(element_number_map[element.first].number);
                input_str += " "
            }
            input_str += "\n"
            input_str += "\n"
            if (this->coordtype == "cartesian") {
                input_str += "xangst" + std::to_string(n) + "\n";
                for (auto atom : this->xyz.atoms) {
                    input_str += std::to_string(atom.x) + " " + std::to_string(atom.y) + " " + std::to_string(atom.z);
                }
            } else if (self.coordtype == "reduced") {
               arma::mat latcell(3, 3);
               latcell.row(0) = arma::conv_to(arma::rowvec)::from(this->xyz.cell[0]);
               latcell.row(1) = arma::conv_to(arma::rowvec)::from(this->xyz.cell[1]);
               latcell.row(2) = arma::conv_to(arma::rowvec)::from(this->xyz.cell[2]);
               auto convmat = arma::inv(arma::trans(latcell));
               arma::mat crystal_coord(this->xyz.natom(), 3);
               arma::vec tmp_vec(3);
               for (int i = 0; i < this->xyz.natom(); i++) {
                   tmp_vec.at(0) = this->xyz.atoms[i].x;
                   tmp_vec.at(1) = this->xyz.atoms[i].y;
                   tmp_vec.at(2) = this->xyz.atoms[i].z;
                   crystal_coord.row(i) = convmat * tmp_vec;
               }
               input_str += "xred" + std::to_string(n) + "\n";
               for (int k = 0; k < crystal_coord.n_rows; k++) {
                   input_str += std::to_string(crystal_coord.at(k, 0)) + " " + std::to_string(crystal_coord.at(k, 1)) + " " + std::to_string(crystal_coord.at(k, 2)) + "\n";
               }
            }
            input_str += "\n";
        } else {
            input_str += "acell" + " 1 1 1 angstrom\n\n";
            input_str += "rprim" + "\n";
            for (int i = 0; i < 3; i++) {
                input_str += std::to_string(cell[i][0]) + " " + std::to_string(cell[i][1]) + " " + std::to_string(cell[i][2]) +"\n";
            }
            input_str += "\n";
            input_str += "ntypat" + " " + std::to_string(this->xyz.nspecies) + "\n\n";
            input_str += "natom" + " " + std::to_string(this->xyz.natom()) + "\n\n";

            input_str += "typat" + " " + "\n";
            // abinit 不允许输入文件列数超过264, 因此如果原子数太多
            // 这里的typat要分多行列出
            // 利用余数设置如果一行超过30个原子就换行
            int i = 0;
            for (auto atom : this->xyz.atoms) {
                input_str += std::to_string(this->xyz.specie_labels[atom.name]) + " ";
                if (i % 30 == 29) {
                    input_str += "\n";
                }
                i += 1;
            }
            input_str += "\n\n";

            input_str += "znucl" + "\n";

            for (auto element :this->xyz.specie_labels) {
                input_str += std::to_string(element_number_map[element.first].number);
                input_str += " "
            }
            input_str += "\n"
            input_str += "\n"
            if (this->coordtype == "cartesian") {
                input_str += "xangst" + std::to_string(n) + "\n";
                for (auto atom : this->xyz.atoms) {
                    input_str += std::to_string(atom.x) + " " + std::to_string(atom.y) + " " + std::to_string(atom.z);
                }
            } else if (self.coordtype == "reduced") {
               arma::mat latcell(3, 3);
               latcell.row(0) = arma::conv_to(arma::rowvec)::from(this->xyz.cell[0]);
               latcell.row(1) = arma::conv_to(arma::rowvec)::from(this->xyz.cell[1]);
               latcell.row(2) = arma::conv_to(arma::rowvec)::from(this->xyz.cell[2]);
               auto convmat = arma::inv(arma::trans(latcell));
               arma::mat crystal_coord(this->xyz.natom(), 3);
               arma::vec tmp_vec(3);
               for (int i = 0; i < this->xyz.natom(); i++) {
                   tmp_vec.at(0) = this->xyz.atoms[i].x;
                   tmp_vec.at(1) = this->xyz.atoms[i].y;
                   tmp_vec.at(2) = this->xyz.atoms[i].z;
                   crystal_coord.row(i) = convmat * tmp_vec;
               }
               input_str += "xred" + "\n";
               for (int k = 0; k < crystal_coord.n_rows; k++) {
                   input_str += std::to_string(crystal_coord.at(k, 0)) + " " + std::to_string(crystal_coord.at(k, 1)) + " " + std::to_string(crystal_coord.at(k, 2)) + "\n";
               }
            }
            input_str += "\n";

        }

        return input_str
    };
    //
    Crystal* crystal;
    bool status;
    std::string coordtype;
};

} // namespace atomsciflow

#endif // atomsciflow_INCLUDE_ASKIT_ABINIT_SYSTEM_H_
