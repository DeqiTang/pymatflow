/************************************************************************
    > File Name: cp2k.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sun 31 Jan 2021 07:28:30 PM CST
************************************************************************/

#include "atomsciflow/cp2k/cp2k.h"

#include <iostream>

namespace atomsciflow {

Cp2kSectionV1& Cp2k::set_subsys(Crystal crystal) {
    this->crystal = crystal;
    return this->set_subsys();
}

Cp2kSectionV1& Cp2k::set_subsys() {
    auto& subsys = this->sections["force_eval"].add_subsection("subsys"); 
    auto& cell = subsys.add_subsection("cell");
    cell.set_param("a", this->crystal.cell[0]);
    cell.set_param("b", this->crystal.cell[1]);
    cell.set_param("c", this->crystal.cell[2]);
    
    // std::cout << "cell set finished!\n"; 
    auto& coord = subsys.add_subsection("coord");
    std::vector<std::vector<std::string>> matrix_str;
    for (auto atom : this->crystal.atoms) {
        matrix_str.push_back(std::vector<std::string>{
                atom.name,
                std::to_string(atom.x),
                std::to_string(atom.y),
                std::to_string(atom.z),
        });
    }
    //
    coord.section_var.set("", matrix_str);
    //std::cout << "subsys set finished!\n"; 
    return subsys;
}


} // namespace atomsciflow

