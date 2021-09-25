/************************************************************************
    > File Name: abinit.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sun 31 Jan 2021 12:31:04 AM CST
************************************************************************/

#ifndef atomsciflow_INCLUDE_ASKIT_ABINIT_ABINIT_H_
#define atomsciflow_INCLUDE_ASKIT_ABINIT_ABINIT_H_

#include "atomsciflow/abinit/electrons.h"


namespace atomsciflow {

class Abinit {

    public:
    
    Abinit() {};
    ~Abinit() {};
   
    std::string to_string() {
        return this->electrons.to_string(0);
    }


    AbinitElectrons electrons;

    private:

};


} // namespace atomsciflow

#endif // atomsciflow_INCLUDE_ASKIT_ABINIT_ABINIT_H_

