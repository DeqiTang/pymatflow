/************************************************************************
    > File Name: gen_section_v1.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sun 31 Jan 2021 08:26:46 PM CST
************************************************************************/
#ifndef atomsciflow_INCLUDE_ASKIT_CP2K_GEN_SECTION_V1_H_
#define atomsciflow_INCLUDE_ASKIT_CP2K_GEN_SECTION_V1_H_

#include "atomsciflow/cp2k/section_v1.h"

namespace atomsciflow {

namespace cp2k {

namespace gen {

    Cp2kSectionV1 global();

    Cp2kSectionV1 force_eval();

    Cp2kSectionV1 dft();
    

} //namespace gen

} //namespace cp2k

} //namespace atomsciflow


#endif // atomsciflow_INCLUDE_ASKIT_CP2K_GEN_SECTION_V1_H_
