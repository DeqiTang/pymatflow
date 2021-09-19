/************************************************************************
    > File Name: group.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sat 30 Jan 2021 08:29:12 PM CST
************************************************************************/

#include "askit/abinit/group_v2.h"

namespace askit {

using namespace utils;


void AbinitVariableGroupV2::set_param(std::string key, int value) {
    this->remove(key);
    this->params_int[key] = AbinitVariableV2<int>{key, value};
}

void AbinitVariableGroupV2::set_param(std::string key, double value) {

    this->remove(key);
    this->params_double[key] = AbinitVariableV2<double>{key, value};
}

void AbinitVariableGroupV2::set_param(std::string key, std::vector<int> value) {
    this->remove(key);

    this->params_int[key] = AbinitVariableV2<int>{key, value};
}

void AbinitVariableGroupV2::set_param(std::string key, std::vector<double> value) {
    this->remove(key);

    this->params_double[key] = AbinitVariableV2<double>{key, value};
}

    bool AbinitVariableGroupV2::contains(std::string key) {
    if (this->params_int.find(key) == this->params_int.end() && this->params_double.find(key) == this->params_double.end()) {
        return false;
    } else {
        return true;
    }
}

void AbinitVariableGroupV2::set_status(std::string key, bool status) {
    //this->params[key].status = status;
    if (this->contains(key) == false) {
        return;
    } else {
        if (this->params_int.find(key) == this->params_int.end()) {
            this->params_double[key].status = status;
        } else {
            this->params_int[key].status = status;
        }
    }
}

std::string AbinitVariableGroupV2::to_string() {
    return this->to_string(this->n);
}

std::string AbinitVariableGroupV2::to_string(int n) {
    std::string out = "";
    for (auto item : this->params_int) {
        if (this->params_int[item.first].status == false) {
            continue;
        }
        out += this->params_int[item.first].to_string(n) + "\n";
    }
    for (auto item : this->params_double) {
        if (this->params_double[item.first].status = false) {
            continue;
        }
        out += this->params_double[item.first].to_string(n) + "\n";
    }
    return out;
}

void AbinitVariableGroupV2::remove(std::string key) {
    for (auto it = this->params_int.begin(); it != this->params_int.end(); ++it) {
        if (it->first == key) {
            this->params_int.erase(it);
            //return;
        }
    }
    for (auto it = this->params_double.begin(); it != this->params_double.end(); ++it) {
        if (it->first == key) {
            this->params_double.erase(it);
            //return;
        }
    }
}


void AbinitVariableGroupV2::clear() {
    this->params_int.clear();
    this->params_double.clear();
}

template<typename U>
U AbinitVariableGroupV2::get(std::string key) {
    U out;
    for (auto it = this->params_int.begin(); it != this->params_int.end(); ++it) {
        if (it->first == key) {
            //this->params_int[key].to(out);
            out = this->params_int[key].as<U>();
            //return out;
        }
    }
    /*
    for (auto it = this->params_double.begin(); it != this->params_double.end(); ++it) {
        if (it->first == key) {
            this->params_double[key].to(out);
            this->params_double[key].as<U>();
            //return out;
        }   
    }
    */
    return out;
}

template int AbinitVariableGroupV2::get(std::string);
//template double AbinitVariableGroupV2::get(std::string);
template std::vector<int> AbinitVariableGroupV2::get(std::string);

} // namespace askit
