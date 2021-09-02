/************************************************************************
    > File Name: group.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sat 30 Jan 2021 08:29:12 PM CST
************************************************************************/

#include "askit/abinit/group_v1.h"

namespace askit {

    using namespace utils;


    void AbinitVariableGroupV1::set_param(std::string key, int value) {
        this->remove(key);
        this->params[key] = AbinitVariableV1{key, value};
    }

    void AbinitVariableGroupV1::set_param(std::string key, double value) {
        this->remove(key);

        this->params[key] = AbinitVariableV1{key, value};
    }

    void AbinitVariableGroupV1::set_param(std::string key, std::string value) {
        this->remove(key);
        this->params[key] = AbinitVariableV1{key, value};
    }

    void AbinitVariableGroupV1::set_param(std::string key, std::vector<int> value) {
        this->remove(key);

        this->params[key] = AbinitVariableV1{key, value};
    }

    void AbinitVariableGroupV1::set_param(std::string key, std::vector<double> value) {
        this->remove(key);

        this->params[key] = AbinitVariableV1{key, value};
    }

    void AbinitVariableGroupV1::set_param(std::string key, std::vector<std::string> value) {
        this->remove(key);

        this->params[key] = AbinitVariableV1{key, value};
    }

    void AbinitVariableGroupV1::set_param(std::string key, std::vector<std::vector<int>> value) {
        this->remove(key);

        this->params[key] = AbinitVariableV1{key, value};
    }

    void AbinitVariableGroupV1::set_param(std::string key, std::vector<std::vector<double>> value) {
        this->remove(key);

        this->params[key] = AbinitVariableV1{key, value};
    }

    void AbinitVariableGroupV1::set_param(std::string key, std::vector<std::vector<std::string>> value) {
        this->remove(key);

        this->params[key] = AbinitVariableV1{key, value};
    }


    bool AbinitVariableGroupV1::contains(std::string key) {
        if (this->params.find(key) == this->params.end()) {
            return false;
        } else {
            return true;
        }
    }

    void AbinitVariableGroupV1::set_status(std::string key, bool status) {
        if (this->contains(key) == false) {
            return;
        }
        this->params[key].status = status;
    }


    std::string AbinitVariableGroupV1::to_string() {
        return this->to_string(this->n);
    }

    std::string AbinitVariableGroupV1::to_string(int n) {
        std::string out = "";
        for (auto item : this->params) {
            if (this->params[item.first].status == false) {
                continue;
            }
            out += this->params[item.first].to_string(n) + "\n";
        }
        return out;
    }


    void AbinitVariableGroupV1::remove(std::string key) {
        for (auto it = this->params.begin(); it != this->params.end(); ++it) {
            if (it->first == key) {
                this->params.erase(it);
                //return;
            }
        }
    }

    void AbinitVariableGroupV1::clear() {
        this->params.clear();
    }


    template<typename U>
    U AbinitVariableGroupV1::get(std::string key) {
        U out;
        //this->params[key].to(out);
        out = this->params[key].as<U>();
        return out;
    }


    template int AbinitVariableGroupV1::get<int>(std::string);
    template double AbinitVariableGroupV1::get<double>(std::string);
    template std::string AbinitVariableGroupV1::get<std::string>(std::string);
    template std::vector<int> AbinitVariableGroupV1::get<std::vector<int>>(std::string);
    template std::vector<double> AbinitVariableGroupV1::get<std::vector<double>>(std::string);
    template std::vector<std::string> AbinitVariableGroupV1::get<std::vector<std::string>>(std::string);
    template std::vector<std::vector<int>> AbinitVariableGroupV1::get<std::vector<std::vector<int>>>(std::string);
    template std::vector<std::vector<double>> AbinitVariableGroupV1::get<std::vector<std::vector<double>>>(std::string);
    template std::vector<std::vector<std::string>> AbinitVariableGroupV1::get<std::vector<std::vector<std::string>>>(std::string);


} // namespace askit
