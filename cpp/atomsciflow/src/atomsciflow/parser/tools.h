/************************************************************************
    > File Name: tools.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 23 Feb 2021 09:35:00 PM CST
************************************************************************/

#ifndef ATOMSCIFLOW_PARSER_TOOLS_H_
#define ATOMSCIFLOW_PARSER_TOOLS_H_

#include "atomsciflow/base/crystal.h"

namespace atomsciflow {


Crystal read_structure_file(std::string filepath);


} // namesapce atomsciflow

#endif // ATOMSCIFLOW_PARSER_TOOLS_H_
