/************************************************************************
    > File Name: src/parser/tools.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 23 Feb 2021 09:38:15 PM CST
************************************************************************/

#include "atomsciflow/parser/tools.h"

#include <boost/filesystem.hpp>
#include "atomsciflow/parser/cif.h"
#include "atomsciflow/parser/xyz.h"
#include "atomsciflow/base/crystal.h"
#include "atomsciflow/utils.h"


namespace atomsciflow {
    
    namespace filesys = boost::filesystem;

    Crystal read_structure_file(std::string filepath) {
        // read structure file
        Crystal crystal;
        filesys::path in_path(filepath);
        //std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;

        if (in_path.has_extension() && in_path.extension().string() == ".cif") {
            atomsciflow::read_cif_file(&crystal, filepath);
        } else if (in_path.has_extension() && in_path.extension().string() == ".xyz") {
            atomsciflow::read_xyz_file(&crystal, filepath);
        }
        
        return crystal;
    }

} // namespace atomsciflow
