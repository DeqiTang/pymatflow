/************************************************************************
    > File Name: section_v1.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sun 31 Jan 2021 07:46:25 PM CST
************************************************************************/

#ifndef ATOMSCIFLOW_CP2K_SECTION_V1_H_
#define ATOMSCIFLOW_CP2K_SECTION_V1_H_

#include <map>

#include "atomsciflow/abinit/variable_v1.h"


namespace atomsciflow {

    using Cp2kVariableV1 = AbinitVariableV1;


class Cp2kSectionV1 {
    public:
    Cp2kSectionV1() {};
    explicit Cp2kSectionV1(std::string name) { this->name = name; }
    ~Cp2kSectionV1() {};


    std::string to_string();
    std::string to_string(std::string indent);

    //void add_subsection(std::string);
    //void add_subsection(std::string, Cp2kSectionV1);
    Cp2kSectionV1& add_subsection(std::string);
    Cp2kSectionV1& add_subsection(std::string, Cp2kSectionV1);

    void remove_subsection(std::string);

    void set_param(std::string key, int value);
    void set_param(std::string key, double value);
    void set_param(std::string key, std::string value);
    void set_param(std::string key, std::vector<int> value);
    void set_param(std::string key, std::vector<double> value);
    void set_param(std::string key, std::vector<std::string> value);

    void set_param(std::string key, std::vector<std::vector<int>> value);
    void set_param(std::string key, std::vector<std::vector<double>> value);
    void set_param(std::string key, std::vector<std::vector<std::string>> value);

    bool contains(std::string key);
    void set_status(std::string key, bool status);
    void remove(std::string key);

    void clear();

    template<typename U>
    U get(std::string key);

    std::string name = "unknown";
    std::string section_parameter;
    Cp2kVariableV1 section_var;
    
    bool status = true;

    private:
    std::map<std::string, Cp2kVariableV1> params;
    std::map<std::string, Cp2kSectionV1> subsections;
};


} // namespace atomsciflow



#endif // ATOMSCIFLOW_CP2K_SECTION_V1_H_

