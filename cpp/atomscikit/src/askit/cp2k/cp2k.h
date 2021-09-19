/************************************************************************
    > File Name: cp2k.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sun 31 Jan 2021 07:24:18 PM CST
************************************************************************/
#ifndef ATOMSCIKIT_INCLUDE_ASKIT_CP2K_CP2K_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_CP2K_CP2K_H_

#include <string>

#include "askit/cp2k/gen_section_v1.h"
#include "askit/base/crystal.h"

namespace askit {

class Cp2k {

public:

    Cp2k() {
        this->sections["global"] = cp2k::gen::global();
        this->sections["force_eval"] = cp2k::gen::force_eval();
    };
    ~Cp2k() {};

    std::string to_string() {
        std::string out = "";
        
        for (auto item : this->sections) {
            out += this->sections[item.first].to_string("  ") + "\n";
            out += "\n";
        }
        return out;
    }


    Cp2kSectionV1& set_subsys(Crystal crystal);

    std::map<std::string, Cp2kSectionV1> sections;
    Crystal crystal;

private:
    Cp2kSectionV1& set_subsys();

};



} // namespace askit


#endif // ATOMSCIKIT_INCLUDE_ASKIT_CP2K_CP2K_H_
