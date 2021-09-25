/************************************************************************
    > File Name: qe_namelist_v1.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 02 Feb 2021 05:10:31 PM CST
************************************************************************/

#include "atomsciflow/qe/namelist_v1.h"

namespace atomsciflow {

    void QeNamelistV1::set_param(std::string key, int value) {
        this->remove(key);
        this->params[key] = QeVariableV1{key, value};
    }

    void QeNamelistV1::set_param(std::string key, double value) {
        this->remove(key);

        this->params[key] = QeVariableV1{key, value};
    }

    void QeNamelistV1::set_param(std::string key, std::string value) {
        this->remove(key);
        this->params[key] = QeVariableV1{key, value};
    }

    void QeNamelistV1::set_param(std::string key, std::vector<int> value) {
        this->remove(key);

        this->params[key] = QeVariableV1{key, value};
    }

    void QeNamelistV1::set_param(std::string key, std::vector<double> value) {
        this->remove(key);

        this->params[key] = QeVariableV1{key, value};
    }

    void QeNamelistV1::set_param(std::string key, std::vector<std::string> value) {
        this->remove(key);

        this->params[key] = QeVariableV1{key, value};
    }

    void QeNamelistV1::set_param(std::string key, std::vector<std::vector<int>> value) {
        this->remove(key);

        this->params[key] = QeVariableV1{key, value};
    }

    void QeNamelistV1::set_param(std::string key, std::vector<std::vector<double>> value) {
        this->remove(key);

        this->params[key] = QeVariableV1{key, value};
    }

    void QeNamelistV1::set_param(std::string key, std::vector<std::vector<std::string>> value) {
        this->remove(key);

        this->params[key] = QeVariableV1{key, value};
    }


    bool QeNamelistV1::contains(std::string key) {
        if (this->params.find(key) == this->params.end()) {
            return false;
        } else {
            return true;
        }
    }

    void QeNamelistV1::set_status(std::string key, bool status) {
        if (this->contains(key) == false) {
            return;
        }
        this->params[key].status = status;
    }


    void QeNamelistV1::remove(std::string key) {
        for (auto it = this->params.begin(); it != this->params.end(); ++it) {
            if (it->first == key) {
                this->params.erase(it);
                //return;
            }
        }
    }

    void QeNamelistV1::clear() {
        this->params.clear();
    }

    std::string QeNamelistV1::to_string() {
        return this->to_string(" ");
    }

    std::string QeNamelistV1::to_string(std::string indent) {
        std::string out = "";
        if (0 == this->type) {
            out += "&" + this->name + "\n";
        } else if (1 == this->type) {
            out += this->name + "\n";
        } else {
            out += "! Warning: unknow namelist type, take as 0 by default\n";
            out += "&" + this->name + "\n";
        }
        for (auto item : this->params) {
            if (false == this->params[item.first].status) {
                continue;
            }
            out += indent + this->params[item.first].to_string() + "\n";
        }

        if (0 == this->type) {
            out += "/\n";
        } else if (1 == this->type) {
            out += "\n";
        } else {
            out += "! Warning: unknow namelist type, take as 0 by default\n";
            out += "/\n";
        }

        return out;
    }


    template<typename U>
    U QeNamelistV1::get(std::string key) {
        U out;
        //this->params[key].to(out);
        out = this->params[key].as<U>();
        return out;
    }


    template int QeNamelistV1::get<int>(std::string);
    template double QeNamelistV1::get<double>(std::string);
    template std::string QeNamelistV1::get<std::string>(std::string);
    template std::vector<int> QeNamelistV1::get<std::vector<int>>(std::string);
    template std::vector<double> QeNamelistV1::get<std::vector<double>>(std::string);
    template std::vector<std::string> QeNamelistV1::get<std::vector<std::string>>(std::string);
    template std::vector<std::vector<int>> QeNamelistV1::get<std::vector<std::vector<int>>>(std::string);
    template std::vector<std::vector<double>> QeNamelistV1::get<std::vector<std::vector<double>>>(std::string);
    template std::vector<std::vector<std::string>> QeNamelistV1::get<std::vector<std::vector<std::string>>>(std::string);


} //namespace atomsciflow
