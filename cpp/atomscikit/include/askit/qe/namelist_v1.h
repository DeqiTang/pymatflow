/************************************************************************
    > File Name: namelist_v1.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 02 Feb 2021 04:53:09 PM CST
************************************************************************/

#ifndef ATOMSCIKIT_INCLUDE_ASKIT_QE_NAMELIST_V1_H_
#define ATOMSCIKIT_INCLUDE_ASKIT_QE_NAMELIST_V1_H_

#include <map>

#include "askit/abinit/variable_v1.h"

namespace askit {

using QeVariableV1 = AbinitVariableV1; 

class QeNamelistV1 {
    public:
    QeNamelistV1() {};
    explicit QeNamelistV1(std::string name) { this->name = name; }
    ~QeNamelistV1() {};

    void set_param(std::string key, int value);
    void set_param(std::string key, double value);
    void set_param(std::string key, std::string value);
    void set_param(std::string key, std::vector<int> value);
    void set_param(std::string key, std::vector<double> value);
    void set_param(std::string key, std::vector<std::string> value);

    void set_param(std::string key, std::vector<std::vector<int>> value);
    void set_param(std::string key, std::vector<std::vector<double>> value);
    void set_param(std::string key, std::vector<std::vector<std::string>> value);

    std::string to_string();
    std::string to_string(std::string indent);

    bool contains(std::string key);
    void set_status(std::string key, bool status);
    void remove(std::string key);

    void clear();

    template<typename U>
    U get(std::string key);

    std::string name = "unknown";
    bool status = true;

    /*
     * type:
     *  0 (with & and /):
     *      &name
     *      /
     *  1 (without & and /):
     *      name
     *
     */
    int type = 0;

    private:
    std::map<std::string, QeVariableV1> params; 
}; 


} // namespace askit

#endif // ATOMSCIKIT_INCLUDE_ASKIT_QE_NAMELIST_V1_H_
