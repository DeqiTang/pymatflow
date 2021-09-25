/************************************************************************
    > File Name: utils.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Mon 25 Jan 2021 08:22:31 PM CST
************************************************************************/

#ifndef atomsciflow_INCLUDE_ASKIT_ABINIT_UTILS_H_
#define atomsciflow_INCLUDE_ASKIT_ABINIT_UTILS_H_

#include <string>
#include <vector>

namespace atomsciflow {

namespace utils {

// return std::to_string(n) if n > 0, else return ""
inline std::string n_to_string(int n) {
    if (n > 0) {
        return std::to_string(n);
    } else {
        return "";
    }
}


} // namespace utils
} // namespace atomsciflow

#endif // atomsciflow_INCLUDE_ASKIT_ABINIT_UTILS_H_
