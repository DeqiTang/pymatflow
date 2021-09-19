/************************************************************************
    > File Name: electrons.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Mon 25 Jan 2021 07:24:49 PM CST
************************************************************************/

#ifndef ATOMSCIKIT_INCLUDE_ASKIT_ABINIT_ELECTRONS_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_ABINIT_ELECTRONS_H_


#include <string>
#include <map>


//#include "askit/abinit/utils.h"
#include "askit/base/kpath.h"
#include "askit/abinit/group_v1.h"
//#include "askit/abinit/group_v2.h"

/*
 * in control of electrons step related parameters
 */


namespace askit {

class Kpoints {
    /*
    default behavior:
        use koptopt = 1 and ngkpt = [1, 1, 1]
    Note:
        I found the setting of kpoints in abinit is versatile
        and I haven't grasped the command of most of them.
        So at present only use the simplest way namely using
        kptopt == 1.
    Reference:
        https://docs.abinit.org/variables/basic/#kptopt
        https://docs.abinit.org/variables/gstate/#kptbounds
    */
    public:

    Kpoints() {
        this->vars = new AbinitVariableGroupV1();
        //this->vars = new AbinitVariableGroupV2;

        this->incharge.push_back("kptopt");
        this->incharge.push_back("ngkpt");
        this->incharge.push_back("nshiftk");
        this->incharge.push_back("shiftk");
        this->basic_setting();
    }

    ~Kpoints() {
        delete this->vars;
    }


    void set_param(std::string key, int value);

    void set_param(std::string key, double value);

    void set_param(std::string key, std::vector<int> value); 

    void set_param(std::string key, std::vector<double> value); 

    void basic_setting() {
        this->set_param("kptopt", 1);
        this->set_param("ngkpt", std::vector<int>{1, 1, 1});
    }

    std::string to_string(int n);

    void set_band(Kpath kpath) {
        this->kpath = kpath;
        this->set_param("kptopt", -(this->kpath.nkpoint - 1));
        //
    }


    std::vector<std::string> incharge;

    //AbinitVariableGroupBase* vars;
    AbinitVariableGroupV1* vars;
    //AbinitVariableGroupV2* vars;
    //
    Kpath kpath;
        
};





class AbinitElectrons {
    /*
     */
    public:
    
    AbinitElectrons() {
        this->vars = new AbinitVariableGroupV1();       
        //this->vars = new AbinitVariableGroupV2;

        this->incharge.push_back("ecut");
        this->incharge.push_back("ixc");
        this->incharge.push_back("nstep");
        this->incharge.push_back("diemac");
        this->incharge.push_back("iscf");
        this->incharge.push_back("toldfe");
        this->incharge.push_back("tolwfr");
        this->incharge.push_back("toldff");
        this->incharge.push_back("tolrff");
        this->incharge.push_back("tolvrs");
        this->incharge.push_back("occopt");
        this->incharge.push_back("nband");
        this->incharge.push_back("occ");
        this->incharge.push_back("wtk");
        this->incharge.push_back("prtden");
        this->incharge.push_back("prtdos");

        this->status = true;

        this->basic_setting();
    }

    ~AbinitElectrons() {
        delete this->vars;
    }

    std::string to_string(int n);

    void check_all_params();

    void use_tol(std::string tol, double value);

    void check_scf_criteria();

    void check_kpoints();
    
    void dft_plus_u();

    void basic_setting();
    

    void set_scf_nscf(std::string mode);

    std::vector<std::string> incharge;
    bool status;
    Kpoints kpoints;
    
    //AbinitVariableGroupBase* vars;
    AbinitVariableGroupV1* vars;
    //AbinitVariableGroupV2* vars;
};





}// namespace askit

#endif // ATOMSCIKIT_INCLUDE_ASKIT_ABINIT_ELECTRONS_H_

