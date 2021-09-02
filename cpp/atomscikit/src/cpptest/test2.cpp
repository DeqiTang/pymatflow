#include "cpptest/test2.h"

#include <string>
#include <iostream>

std::string log() {
    std::cout << "test print" << std::endl;
    return std::string("test log string");
}
