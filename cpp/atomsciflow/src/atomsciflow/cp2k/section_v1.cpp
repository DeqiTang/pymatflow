/************************************************************************
    > File Name: section_v1.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sun 31 Jan 2021 07:57:17 PM CST
************************************************************************/
#include "atomsciflow/cp2k/section_v1.h"

namespace atomsciflow {

std::string Cp2kSectionV1::to_string() {

    if (this->status == false) {
        return "";
    }

    std::string out = "";
    
    out += "&" + this->name + " " + this->section_parameter + "\n";

    out += this->section_var.to_string() + "\n";

    for (auto item : this->params) {
        out += this->params[item.first].to_string() + "\n";
    }

    for (auto item : this->subsections) {
        out += "\n";
        out += this->subsections[item.first].to_string();
        out += "\n";
    }

    out += "&end " + this->name;
    return out;
}

std::string Cp2kSectionV1::to_string(std::string indent) {

    if (false == this->status) {
        return "";
    }

    std::string out = "";
    
    out += indent + "&" + this->name + " " + this->section_parameter + "\n";

    out += this->section_var.to_string("same-line", indent + indent) + "\n";

    for (auto item : this->params) {
        out += indent + indent + this->params[item.first].to_string("same-line", "") + "\n";
    }

    for (auto item : this->subsections) {
        //out += "\n";
        out += this->subsections[item.first].to_string(indent + indent);
        out += "\n";
    }
   
    out += indent + "&end " + this->name;
    
    return out;
}

Cp2kSectionV1& Cp2kSectionV1::add_subsection(std::string key) {
    //this->subsections[key] = Cp2kSectionV1{};
    //this->subsections[key].name = key;
    this->subsections[key] = Cp2kSectionV1{key};
    return this->subsections[key];
}


Cp2kSectionV1& Cp2kSectionV1::add_subsection(std::string key, Cp2kSectionV1 section) {
    this->subsections[key] = section;
    this->subsections[key].name = key;
    return this->subsections[key];
}


void Cp2kSectionV1::remove_subsection(std::string key) {
    for (auto it = this->subsections.begin(); it != this->subsections.end(); ++it) {
        if (it->first == key) {
            this->subsections.erase(it);
        }
    }
}


void Cp2kSectionV1::set_param(std::string key, int value) {
    this->remove(key);
    this->params[key] = AbinitVariableV1{key, value};
}

void Cp2kSectionV1::set_param(std::string key, double value) {
    this->remove(key);
    this->params[key] = AbinitVariableV1{key, value};
}

void Cp2kSectionV1::set_param(std::string key, std::string value) {
    this->remove(key);
    this->params[key] = AbinitVariableV1{key, value};
}

void Cp2kSectionV1::set_param(std::string key, std::vector<int> value) {
    this->remove(key);

    this->params[key] = AbinitVariableV1{key, value};
}

void Cp2kSectionV1::set_param(std::string key, std::vector<double> value) {
    this->remove(key);

    this->params[key] = AbinitVariableV1{key, value};
}

void Cp2kSectionV1::set_param(std::string key, std::vector<std::string> value) {
    this->remove(key);

    this->params[key] = AbinitVariableV1{key, value};
}

void Cp2kSectionV1::set_param(std::string key, std::vector<std::vector<int>> value) {
    this->remove(key);

    this->params[key] = AbinitVariableV1{key, value};
}

void Cp2kSectionV1::set_param(std::string key, std::vector<std::vector<double>> value) {
    this->remove(key);

    this->params[key] = AbinitVariableV1{key, value};
}

void Cp2kSectionV1::set_param(std::string key, std::vector<std::vector<std::string>> value) {
    this->remove(key);

    this->params[key] = AbinitVariableV1{key, value};
}

bool Cp2kSectionV1::contains(std::string key) {
    if (this->params.find(key) == this->params.end()) {
        return false;
    } else {
        return true;
    }
}

void Cp2kSectionV1::set_status(std::string key, bool status) {
    if (this->contains(key) == false) {
        return;
    }
    this->params[key].status = status;
}

void Cp2kSectionV1::remove(std::string key) {
    for (auto it = this->params.begin(); it != this->params.end(); ++it) {
        if (it->first == key) {
            this->params.erase(it);
            //return;
        }
    }
}

void Cp2kSectionV1::clear() {
    this->params.clear();
}


template<typename U>
U Cp2kSectionV1::get(std::string key) {
    U out;
    //this->params[key].to(out);
    out = this->params[key].as<U>();
    return out;
}


template int Cp2kSectionV1::get<int>(std::string);
template double Cp2kSectionV1::get<double>(std::string);
template std::string Cp2kSectionV1::get<std::string>(std::string);
template std::vector<int> Cp2kSectionV1::get<std::vector<int>>(std::string);
template std::vector<double> Cp2kSectionV1::get<std::vector<double>>(std::string);
template std::vector<std::string> Cp2kSectionV1::get<std::vector<std::string>>(std::string);
template std::vector<std::vector<int>> Cp2kSectionV1::get<std::vector<std::vector<int>>>(std::string);
template std::vector<std::vector<double>> Cp2kSectionV1::get<std::vector<std::vector<double>>>(std::string);
template std::vector<std::vector<std::string>> Cp2kSectionV1::get<std::vector<std::vector<std::string>>>(std::string);




} //namespace atomsciflow
