/************************************************************************
    > File Name: gen_force_eval.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sun 31 Jan 2021 08:47:08 PM CST
************************************************************************/

#include "atomsciflow/cp2k/gen_section_v1.h"

#include <iostream>

namespace atomsciflow {

namespace cp2k {

namespace gen {


Cp2kSectionV1 dft() {
    Cp2kSectionV1 out{"dft"};
    
    out.set_param("basis_set_file_name", "basis_set");
    out.set_param("potential_file_name", "gth_poetntials");

    auto& qs = out.add_subsection("qs");
    qs.set_param("eps_default", "1.0e-14");
    
    auto& mgrid = out.add_subsection("mgrid");
    mgrid.set_param("ngrids", 4);
    mgrid.set_param("cutoff", 100);
    mgrid.set_param("rel_cutoff", 60);
    

    auto& xc = out.add_subsection("xc");
    auto& xc_functional = xc.add_subsection("xc_functional");
    xc_functional.section_parameter = "pbe";


    auto& scf = out.add_subsection("scf");
    scf.set_param("scf_guess", "atomic");
    scf.set_param("eps_scf", "1.0e-7");
    scf.set_param("max_scf", 100);
    auto& diag = scf.add_subsection("diagonalization");
    diag.set_param("algorithm", "standard");

    auto& mixing = out.add_subsection("mixing");
    //mixing.value = "t";
    mixing.set_param("method", "broyden_mixing");
    mixing.set_param("alpha", 0.4);
    mixing.set_param("nbroyden", 8);

    return out;
}

Cp2kSectionV1 force_eval() {
    Cp2kSectionV1 out{"force_eval"};

    out.set_param("method", "quickstep");

    out.add_subsection("dft", dft());

    auto& print = out.add_subsection("print");  
    print.add_subsection("forces").section_parameter = "on";

    return out;
}


} // namesapce gen 

} // namespace cp2k

} // namespace atomsciflow
