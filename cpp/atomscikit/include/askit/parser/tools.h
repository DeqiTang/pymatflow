/************************************************************************
    > File Name: tools.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 23 Feb 2021 09:35:00 PM CST
************************************************************************/

#ifndef ATOMSCIKIT_INCLUDE_ASKIT_PARSER_TOOLS_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_PARSER_TOOLS_H_

#include "askit/base/crystal.h"

namespace askit {


Crystal read_structure_file(std::string filepath);


} // namesapce askit

#endif // ATOMSCIKIT_INCLUDE_ASKIT_PARSER_TOOLS_H_
