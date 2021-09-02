/************************************************************************
    > File Name: pw.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 02 Feb 2021 05:56:20 PM CST
************************************************************************/

#include "askit/qe/pw.h"

namespace askit {

std::string QePw::to_string() {

    return this->to_string("");    
}

std::string QePw::to_string(std::string indent) {

    std::string out = "";
    for (auto item : this->namelists) {
        if (false == this->namelists[item.first].status) {
            continue;
        } else {
            out += this->namelists[item.first].to_string(indent);
            out += "\n";
        }
    }
    return out;
}

} // namespcae askit
