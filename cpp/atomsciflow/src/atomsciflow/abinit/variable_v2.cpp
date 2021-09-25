/************************************************************************
    > File Name: variable_v2.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sat 30 Jan 2021 01:11:54 PM CST
************************************************************************/

//#include "atomsciflow/abinit/utils.h"
#include "atomsciflow/abinit/variable_v2.h"

#include <string>
//#include <vector>

namespace atomsciflow {
    
    using namespace utils;


    template<typename T>
    void AbinitVariableV2<T>::set(std::string key, T value) {
        this->key = key;
        this->value.clear();
        this->value.push_back(std::vector<T>{value});
    }

    template<typename T>
    void AbinitVariableV2<T>::set(std::string key, std::vector<T> value) {
        this->key = key;
        this->value.clear();
        this->value.push_back(value);
    }

    template<typename T>
    void AbinitVariableV2<T>::set(std::string key, std::vector<std::vector<T> > value) {
        this->key = key;
        this->value = value;
    }

    /*
     * if we do not specify the object type explicityly
     * we can only export to calss template type T
     * and that's not enough
     */
    template<typename T>
    void AbinitVariableV2<T>::to(T& value) {
        value = this->value[0][0];
    }

    template<typename T>
    void AbinitVariableV2<T>::to(std::vector<T>& value) {
        value = this->value[0];
    }

    template<typename T>
    void AbinitVariableV2<T>::to(std::vector<std::vector<T>>& value) {
        value = this->value;
    }



    template<typename T>
    std::string AbinitVariableV2<T>::to_string() {
        return this->to_string(this->n);
    }


    template<typename T>
    std::string AbinitVariableV2<T>::to_string(int n) {
        std::string out = "";
        if (this->value.size() == 1) {
            if (this->value[0].size() == 1) {
                out += this->key + n_to_string(n) + " " + std::to_string(this->value[0][0]);
            } else {
                out += this->key + n_to_string(n);
                for (int i = 0; i < this->value[0].size(); i ++) {
                    out += " " + std::to_string(this->value[0][i]);
                }
            }
        } else {
            out += this->key + n_to_string(n) + "\n";
            for (int i = 0; i < this->value.size(); i++) {
                for (int j = 0; j < this->value[i].size(); j++) {
                    out += " " + std::to_string(this->value[i][j]);
                }
                out += "\n";
            }
        }
        return out;
    }

    // explicit template instantiation
    template class AbinitVariableV2<int>;
    template class AbinitVariableV2<double>;


} // namespace atomsciflow
