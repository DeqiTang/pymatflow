/************************************************************************
    > File Name: variable_v1.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Mon 25 Jan 2021 08:22:31 PM CST
************************************************************************/

#ifndef atomsciflow_INCLUDE_atomsciflow_ABINIT_VARIABLE_V1_H_
#define atomsciflow_INCLUDE_atomsciflow_ABINIT_VARIABLE_V1_H_

#include <string>
#include <vector>

#include "atomsciflow/abinit/utils.h"

namespace atomsciflow {

class AbinitVariableV1 {
    public:
        AbinitVariableV1() {};
        ~AbinitVariableV1() {};

        AbinitVariableV1(std::string key, int value) {
            this->set(key, value);
        }

        AbinitVariableV1(std::string key, double value) {
            this->set(key, value);
        }

        AbinitVariableV1(std::string key, std::string value) {
            this->set(key, value);
        }

        AbinitVariableV1(std::string key, std::vector<int> value) {
            this->set(key, value);
        }

        AbinitVariableV1(std::string key, std::vector<double> value) {
            this->set(key, value);
        }

        AbinitVariableV1(std::string key, std::vector<std::string> value) {
            this->set(key, value);
        }

        AbinitVariableV1(std::string key, std::vector<std::vector<int> > value) {
            this->set(key, value);
        }

        AbinitVariableV1(std::string key, std::vector<std::vector<double> > value) {
            this->set(key, value);
        }

        AbinitVariableV1(std::string key, std::vector<std::vector<std::string> > value) {
            this->set(key, value);
        }

        void set(std::string key, int value);
        void set(std::string key, double value);
        void set(std::string key, std::string value);
        void set(std::string key, std::vector<int> value);
        void set(std::string key, std::vector<double> value);
        void set(std::string key, std::vector<std::string> value);
        void set(std::string key, std::vector<std::vector<int>> value);
        void set(std::string key, std::vector<std::vector<double>> value);
        void set(std::string key, std::vector<std::vector<std::string>> value);

        // with implicit instantiation
        template<typename U>
        U as() {
            U out;
            this->to(out);
            return out;
        }

        void to(int& value);
        void to(double& value);
        void to(std::string& value);
        void to(std::vector<int>& value);
        void to(std::vector<double>& value);
        void to(std::vector<std::string>& value);
        void to(std::vector<std::vector<int>>& value);
        void to(std::vector<std::vector<double>>& value);
        void to(std::vector<std::vector<std::string>>& value);

        void set_n(int n) { this->n = n; }

        std::string to_string();
        std::string to_string(int n);
        std::string to_string(std::string layout, std::string indent);
        std::string to_string(std::string layout, std::string indent, int n);


        std::string key;
        int n;
        std::vector<std::vector<std::string>> value;

        bool status = true;
    private:
};



} // namespace atomsciflow
#endif // atomsciflow_INCLUDE_atomsciflow_ABINIT_VARIABLE_V1_H_
