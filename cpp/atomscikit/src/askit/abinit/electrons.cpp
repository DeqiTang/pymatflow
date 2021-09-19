/************************************************************************
    > File Name: ../../../src/abinit/electrons.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sat 30 Jan 2021 02:37:52 PM CST
************************************************************************/

#include "askit/abinit/electrons.h"

namespace askit {

void Kpoints::set_param(std::string key, int value) {
    this->vars->set_param(key, value);
}

void Kpoints::set_param(std::string key, double value) {
    this->vars->set_param(key, value);
}

void Kpoints::set_param(std::string key, std::vector<int> value) {
    this->vars->set_param(key, value);
}

void Kpoints::set_param(std::string key, std::vector<double> value) {
    this->vars->set_param(key, value);
}


std::string Kpoints::to_string(int n = 0) {
    /*
    :return input_str is the string of all the set params
    */
    std::string out = "";

    out += this->vars->to_string(n);

    // end
    //
    out += "\n";

    return out;
}


// AbinitElectrons
// ---------------------------------------

std::string AbinitElectrons::to_string(int n = 0) {
    /*
     *
     :return input_str is the string of all the set params
    */
    std::string out = "";

    out += this->vars->to_string(n);

    out += "\n";
    return out;
}

void AbinitElectrons::check_all_params() {
    //
    //
}


void AbinitElectrons::use_tol(std::string tol, double value) {
    //
    std::vector<std::string> tols{"toldfe", "tolwfr", "toldff", "tolrff", "tolvrs"};
    for (auto item : tols) {
        if (item != tol && this->vars->contains(item)) {
            this->vars->set_status(item, false);
        }
    }
    this->vars->set_param(tol, value);
    this->vars->set_status(tol, true);
}




void AbinitElectrons::set_scf_nscf(std::string mode="scf") {
    //
    if (mode == "scf") {
        if (this->vars->contains("usepaw") && this->vars->get<int>("usepaw") == 1) {
            this->vars->set_param("iscf", 17);
        } else {
            this->vars->set_param("iscf", 7);
        }
        this->vars->set_param("prtden", 1);
    }

    if (mode == "nscf") {
        this->vars->set_param("iscf", -3);
        this->vars->set_param("nstep", 0);
    }
}


void AbinitElectrons::check_scf_criteria() {
    /*
     * there are six kinds of criteria to judge convergence of
     * scf in abinit, but one and only one of them can be used
     * in a run. here we make a check
     */
}
   
void AbinitElectrons::check_kpoints() {
    // there should be no kpoints related setting in self.params
}

void AbinitElectrons::dft_plus_u() {
    // 需要使用paw赝势, 在files文件中进行设置
    // 

}


void AbinitElectrons::basic_setting() {
    this->vars->set_param("ecut", 15);
    this->vars->set_param("occopt", 3);
    this->vars->set_param("nstep", 100);
    this->vars->set_param("diemac", 2.0);
    this->vars->set_param("ixc", 11);
    this->use_tol("tolvrs", 1.0e-18);
}


} // namespace askit
