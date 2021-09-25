/************************************************************************
    > File Name: tools.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 23 Feb 2021 09:35:00 PM CST
************************************************************************/

#ifndef atomsciflow_INCLUDE_atomsciflow_PARSER_TOOLS_H_
#define atomsciflow_INCLUDE_atomsciflow_PARSER_TOOLS_H_

#include "atomsciflow/base/crystal.h"

namespace atomsciflow {


Crystal read_structure_file(std::string filepath);


} // namesapce atomsciflow

#endif // atomsciflow_INCLUDE_atomsciflow_PARSER_TOOLS_H_
