/************************************************************************
    > File Name: gen_namelist_v1.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 02 Feb 2021 05:30:54 PM CST
************************************************************************/

#ifndef atomsciflow_INCLUDE_atomsciflow_QE_GEN_NAMELIST_V1_H_
#define atomsciflow_INCLUDE_atomsciflow_QE_GEN_NAMELIST_V1_H_

#include "atomsciflow/qe/namelist_v1.h"

namespace atomsciflow {

namespace qe {
namespace gen {      

QeNamelistV1 control();

QeNamelistV1 system();

QeNamelistV1 electrons();

QeNamelistV1 ions();

QeNamelistV1 cell();


} // namespace gen
} // namespace qe
} // namespace atomsciflow


#endif // atomsciflow_INCLUDE_atomsciflow_QE_GEN_NAMELIST_V1_H_

