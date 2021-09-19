/************************************************************************
    > File Name: gen_global.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sun 31 Jan 2021 08:35:51 PM CST
************************************************************************/

#include "askit/cp2k/gen_section_v1.h" 

namespace askit {

namespace cp2k {

namespace gen {

Cp2kSectionV1 global() {
    Cp2kSectionV1 out{"global"};

    out.set_param("project", "cp2k_job");
    out.set_param("print_level", "low");    
    out.set_param("run_type", "energy_force");
    return out;
}


} // namesapce gen 

} // namespace cp2k

} // namespace askit
