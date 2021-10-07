/************************************************************************
    > File Name: pw.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 02 Feb 2021 05:49:31 PM CST
************************************************************************/


#ifndef ATOMSCIFLOW_QE_PW_H_
#define ATOMSCIFLOW_QE_PW_H_

#include <map>

#include "atomsciflow/qe/gen_namelist_v1.h"

namespace atomsciflow {

class QePw {

    public:
    
    QePw() {
        this->namelists["control"] = qe::gen::control();
        this->namelists["system"] = qe::gen::system();
        this->namelists["electrons"] = qe::gen::electrons();
        this->namelists["ions"] = qe::gen::ions();
        this->namelists["cell"] = qe::gen::cell();
    }

    ~QePw() {};

    std::string to_string();
    std::string to_string(std::string indent);
    
    private:

    std::map<std::string, QeNamelistV1> namelists;
};


} // namespace atomsciflow


#endif // ATOMSCIFLOW_QE_PW_H_
