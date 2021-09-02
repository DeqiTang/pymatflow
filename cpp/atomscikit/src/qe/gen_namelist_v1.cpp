/************************************************************************
    > File Name: gen_namelist_v1.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 02 Feb 2021 05:39:55 PM CST
************************************************************************/

#include "askit/qe/gen_namelist_v1.h"

namespace askit {
 
namespace qe {
namespace gen {

QeNamelistV1 control() {

    QeNamelistV1 out;
    
    out.name = "control";
    out.set_param("calculation", "scf");
    out.set_param("verbosity", "low");

    return out;
}


QeNamelistV1 system() {
    QeNamelistV1 out;

    out.name = "system";
    out.set_param("input_dft", "PBE");

    return out;
}


QeNamelistV1 electrons() {
    QeNamelistV1 out;

    out.name = "electron";
    out.set_param("electron_maxstep", 100);
    out.set_param("conv_thr", "1.0E-6");
    
    return out;
}


QeNamelistV1 ions() {
    QeNamelistV1 out;

    out.name = "ion";
    out.set_param("ion_dynamics", "bfgs");

    return out;
}


QeNamelistV1 cell() {
    QeNamelistV1 out;

    out.name = "cell";
    out.set_param("cell_dynamics", "bfgs");
    return out;
}


} // gen
} // qe
} // askit
