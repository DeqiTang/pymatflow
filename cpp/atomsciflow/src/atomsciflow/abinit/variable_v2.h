/************************************************************************
    > File Name: variable_v2.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Mon 25 Jan 2021 08:22:31 PM CST
************************************************************************/

#ifndef atomsciflow_INCLUDE_atomsciflow_ABINIT_VARIABLE_V2_H_
#define atomsciflow_INCLUDE_atomsciflow_ABINIT_VARIABLE_V2_H_

#include <string>
#include <vector>

#include "atomsciflow/abinit/utils.h"

namespace atomsciflow {


template<typename T>
class AbinitVariableV2 {
    public:
        AbinitVariableV2() {};
        ~AbinitVariableV2() {};

        AbinitVariableV2(std::string key, T value) {
            this->set(key, value);
        }

        AbinitVariableV2(std::string key, std::vector<T> value) {
            this->set(key, value);
        }

        AbinitVariableV2(std::string key, std::vector<std::vector<T> > value) {
            this->set(key, value);
        }

        void set(std::string key, T value);

        void set(std::string key, std::vector<T> value);

        void set(std::string key, std::vector<std::vector<T> > value);

        void set_n(int n) { this->n = n; }

        // with implicit instantiation
        template<typename U>
        U as() {
            U out;
            this->to(out);
            return out;
        }

        /*
         * if we donot specify the object type explicitly
         * we can only export to calass template type T
         * and that's not engough
         */
        void to(T& value);

        void to(std::vector<T>& value);

        void to(std::vector<std::vector<T> >& value);
        
        //void to(int& value);
        //void to(double& value);
        //void to(std::vector<int>& value);
        //void to(std::vector<double>& value);
        //void to(std::vector<std::vector<int>>& value);
        //void to(std::vector<std::vector<double>>& value);
        

        std::string to_string();

        std::string to_string(int n);

        std::string key;
        int n;
        std::vector<std::vector<T> > value;

        bool status = true;
    private:
};

} // namespace atomsciflow

#endif // atomsciflow_INCLUDE_atomsciflow_ABINIT_VARIABLE_V2_H_
