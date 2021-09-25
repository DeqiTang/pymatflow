/************************************************************************
    > File Name: variable_v1.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Sat 30 Jan 2021 01:11:54 PM CST
************************************************************************/

//#include "atomsciflow/abinit/utils.h"
#include "atomsciflow/abinit/variable_v1.h"

#include <string>
//#include <vector>

namespace atomsciflow {

    using namespace utils;
    
    // inline functions
    //
    inline std::string to_string_same_line(const AbinitVariableV1& var, std::string indent, int n) {
        if (false == var.status) {
            return "";
        }
        std::string out = "";
        if (0 == var.value.size()) {
            return out + var.key;
        }
        if (1 == var.value.size()) {
            if (1 == var.value[0].size()) {
                out += indent + var.key + n_to_string(n) + " " + var.value[0][0];
            } else {
                out += indent + var.key + n_to_string(n);
                for (auto item : var.value[0]) {
                    out += " " + item;
                }
            } 
        } else {
            out += indent + var.key + n_to_string(n);
            for (auto val : var.value[0]) {
                out += " " + val;
            }
            out += "\n";
            for (int row = 1; row < var.value.size() - 1; ++row) {
                out += indent;
                for (auto val : var.value[row]) {
                    out += " " + val;
                }
                out += "\n";
            }
            out += indent;
            for (auto val : var.value[var.value.size()-1]) {
                out += " " + val;
            }
        }
        return out;
    }


    inline std::string to_string_second_line(const AbinitVariableV1& var, std::string indent, int n) {
        
        if (false == var.status) {
            return "";
        }
        
        std::string out = "";
        
        if (0 == var.value.size()) {
            return out + var.key;
        }

        if (1 == var.value.size()) {
            if (1 == var.value[0].size()) {
                out += indent + var.key + n_to_string(n) + " " + var.value[0][0];
            } else {
                out += indent + var.key + n_to_string(n) + "\n";

                out += indent;
                for (auto item : var.value[0]) {
                    out += " " + item;
                }
            } 
        } else {
            out += indent + var.key + n_to_string(n) + "\n";
            for (auto row : var.value) {
                out += indent;
                for (auto val : row) {
                    out += " " + val;
                }
                out += "\n";
            }
        }
        return out;
    }


    void AbinitVariableV1::set(std::string key, int value) {
        this->key = key;
        this->value.clear();
        this->value.push_back(std::vector<std::string>{std::to_string(value)});
    }

    void AbinitVariableV1::set(std::string key, double value) {
        this->key = key;
        this->value.clear();
        this->value.push_back(std::vector<std::string>{std::to_string(value)});
    }

    void AbinitVariableV1::set(std::string key, std::string value) {
        this->key = key;
        this->value.clear();
        this->value.push_back(std::vector<std::string>{value});
    }


    void AbinitVariableV1::set(std::string key, std::vector<int> value) {
        this->key = key;
        this->value.clear();
        std::vector<std::string> vec_str;
        for (auto& i : value) {
            vec_str.push_back(std::to_string(i));
        }
        this->value.push_back(vec_str);
    }

    void AbinitVariableV1::set(std::string key, std::vector<double> value) {
        this->key = key;
        this->value.clear();
        std::vector<std::string> vec_str;
        for (auto& i : value) {
            vec_str.push_back(std::to_string(i));
        }
        this->value.push_back(vec_str);

    }

    void AbinitVariableV1::set(std::string key, std::vector<std::string> value) {
        this->key = key;
        this->value.clear();
        std::vector<std::string> vec_str;
        for (auto& val : value) {
            vec_str.push_back(val);
        }
        this->value.push_back(vec_str);

    }

    void AbinitVariableV1::set(std::string key, std::vector<std::vector<int> > value) {
        this->key = key;
        this->value.clear();
        std::vector<std::string> vec_str;
        for (auto& row : value) {
            vec_str.clear();
            for (auto& val : row) {
                vec_str.push_back(std::to_string(val));
            }
            this->value.push_back(vec_str);
        }
    }

    void AbinitVariableV1::set(std::string key, std::vector<std::vector<double> > value) {
        this->key = key;
        this->value.clear();
        std::vector<std::string> vec_str;
        for (auto& row : value) {
            vec_str.clear();
            for (auto& val : row) {
                vec_str.push_back(std::to_string(val));
            }
            this->value.push_back(vec_str);
        }
    }

    void AbinitVariableV1::set(std::string key, std::vector<std::vector<std::string> > value) {
        this->key = key;
        this->value.clear();
        std::vector<std::string> vec_str;
        for (auto& row : value) {
            vec_str.clear();
            for (auto& val : row) {
                vec_str.push_back(val);
            }
            this->value.push_back(vec_str);
        }
    }


    void AbinitVariableV1::to(int& value) {
        value = std::atoi(this->value[0][0].c_str());
    }

    void AbinitVariableV1::to(double& value) {
        value = std::atof(this->value[0][0].c_str());
    }

    void AbinitVariableV1::to(std::string& value) {
        value = this->value[0][0];
    }

    void AbinitVariableV1::to(std::vector<int>& value) {
        value.clear();
        for (auto& val : this->value[0]) {
            value.push_back(std::atoi(val.c_str()));
        }
    }
    
    void AbinitVariableV1::to(std::vector<double>& value) {
        value.clear();
        for (auto& val : this->value[0]) {
            value.push_back(std::atof(val.c_str()));
        }
    }

    void AbinitVariableV1::to(std::vector<std::string>& value) {
        value.clear();
        //for (auto& val : this->value[0]) {
        //    value.push_back(val);
        //}
        value = this->value[0];
    }

    void AbinitVariableV1::to(std::vector<std::vector<int>>& value) {
        value.clear();
        std::vector<int> vec_int;
        for (auto& row : this->value) {
            vec_int.clear();
            for (auto& val : row) {
                vec_int.push_back(std::atoi(val.c_str()));
            }
            value.push_back(vec_int);
        }
    }

    void AbinitVariableV1::to(std::vector<std::vector<double>>& value) {
        value.clear();
        std::vector<double> vec_double;
        for (auto& row : this->value) {
            vec_double.clear();
            for (auto& val : row) {
                vec_double.push_back(std::atof(val.c_str()));
            }
            value.push_back(vec_double);
        }
    }

    void AbinitVariableV1::to(std::vector<std::vector<std::string>>& value) {
        value.clear();
        //std::vector<std::string> vec_str;
        //for (auto& row : this->value) {
        //    vec_double.clear();
        //    for (auto& val : row) {
        //        vec_double.push_back(val);
        //    }
        //    value.push_back(vec_str);
        //}
        value = this->value;
    }


    std::string AbinitVariableV1::to_string() {
        if (false == this->status) {
            return "";
        }
        return this->to_string(0);
    }

    std::string AbinitVariableV1::to_string(int n) {
        if (false == this->status) {
            return "";
        }
        std::string out = "";
        if (9 == this->value.size()) {
            return out + this->key;
        }

        if (this->value.size() == 1) {
            if (this->value[0].size() == 1) {
                out += this->key + n_to_string(n) + " " + this->value[0][0];
            } else {
                out += this->key + n_to_string(n) + "\n";
                for (auto item : this->value[0]) {
                    out += " " + item;
                }
            } 
        } else {
            out += this->key + n_to_string(n) + "\n";
            for (auto row : this->value) {
                for (auto val : row) {
                    out += " " + val;
                }
                out += "\n";
            }
        }
        return out;
    }

    std::string AbinitVariableV1::to_string(std::string layout, std::string indent) {
        if (false == this->status) {
            return "";
        }
        return this->to_string(layout, indent, 0);
    }


    std::string AbinitVariableV1::to_string(std::string layout, std::string indent, int n) {

        /*
         * layout:
         *      "same-line";
         *      "second-line"''
         */
        if (false == this->status) {
            return "";
        }
        if ("same-line" == layout) {
            return to_string_same_line(*this, indent, n);
        } else if ("second-line" == layout) {
            return to_string_second_line(*this, indent, n);
        } else {
            return to_string(n);
        }
    }

    //

} // namespace atomsciflow
