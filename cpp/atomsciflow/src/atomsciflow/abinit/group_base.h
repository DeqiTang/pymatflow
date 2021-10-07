/************************************************************************
    > File Name: group_base.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sat 30 Jan 2021 08:14:26 PM CST
************************************************************************/

#ifndef ATOMSCIFLOW_ABINIT_GROUP_BASE_H_
#define ATOMSCIFLOW_ABINIT_GROUP_BASE_H_

#include <map>
#include <iostream>
#include <iterator>
#include <typeinfo>

#include "atomsciflow/abinit/utils.h"

namespace atomsciflow {

class AbinitVariableGroupBase {

    public:

        AbinitVariableGroupBase() {};
        ~AbinitVariableGroupBase() {};

        virtual void set_param(std::string key, int value);
        virtual void set_param(std::string key, double value);
        virtual void set_param(std::string key, std::vector<int> value);
        virtual void set_param(std::string key, std::vector<double> value);

        virtual std::string to_string();
        virtual std::string to_string(int n);

        virtual bool contains(std::string key);
        virtual void set_status(std::string key, bool status);
        virtual void remove(std::string key);

        int n;
    private:

};


} // namespace atomsciflow

#endif // ATOMSCIFLOW_ABINIT_GROUP_BASE_H_

