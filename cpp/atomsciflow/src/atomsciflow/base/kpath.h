/************************************************************************
    > File Name: kpath.h
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Mon 25 Jan 2021 07:39:31 PM CST
************************************************************************/

#ifndef atomsciflow_INCLUDE_ASKIT_BASE_KPATH_H_
#define atomsciflow_INCLUDE_ASKIT_BASE_KPATH_H_

#include <vector>

/*
 *
 */


namespace atomsciflow {

class Kpath {
    /*
     * the high symmetry k point path used in bands structure calculation
     * in format like this:
     *     [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]
     *     described by a python list, however not in c++
     * if connect_indicator in a kpoint is an integer larger than 0, then it will connect to the following point
     * through the number of kpoints defined by connect_indicator.
     * if connect_indicator in a kpoint is an integer 0, then it will not connect to the following point,
    */
    public:
    Kpath() {};
    ~Kpath() {};

    void add_point(float kx, float ky, float kz, std::string label, int connect_indicator) {
        this->points.push_back(std::vector<double>{kx, ky, kz});
        this->labels.push_back(label);
        this->connect_indicator.push_back(connect_indicator);
        return ;
    }

    std::vector<std::vector<double> > points;
    std::vector<std::string> labels;
    std::vector<int> connect_indicator;
    int nkpoint;
};



} // namespace atomsciflow

#endif // atomsciflow_INCLUDE_ASKIT_BASE_KPATH_H_

