/************************************************************************
    > File Name: group_v1.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sat 30 Jan 2021 08:14:26 PM CST
************************************************************************/

#ifndef ATOMSCIKIT_INCLUDE_ASKIT_ABINIT_GROUP_V1_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_ABINIT_GROUP_V1_H_

#include <map>
#include <iostream>
#include <iterator>
#include <typeinfo>

//#include "askit/abinit/utils.h"
#include "askit/abinit/variable_v1.h"

namespace askit {

class AbinitVariableGroupV1 {//}: public AbinitVariableGroupBase {
    public:
        AbinitVariableGroupV1() {};
        virtual ~AbinitVariableGroupV1() {};

        virtual void set_param(std::string key, int value);
        virtual void set_param(std::string key, double value);
        virtual void set_param(std::string key, std::string value);
        virtual void set_param(std::string key, std::vector<int> value);
        virtual void set_param(std::string key, std::vector<double> value);
        virtual void set_param(std::string key, std::vector<std::string> value);

        virtual void set_param(std::string key, std::vector<std::vector<int>> value);
        virtual void set_param(std::string key, std::vector<std::vector<double>> value);
        virtual void set_param(std::string key, std::vector<std::vector<std::string>> value);

        virtual std::string to_string();
        virtual std::string to_string(int n);

        virtual bool contains(std::string key);
        virtual void set_status(std::string key, bool status);
        virtual void remove(std::string key);

        virtual void clear();

        template<typename U>
        U get(std::string key);

        int n;
    private:
        std::map<std::string, AbinitVariableV1> params;
};


} // namespace askit

#endif // ATOMSCIKIT_INCLUDE_ASKIT_ABINIT_GROUP_V1_H_

