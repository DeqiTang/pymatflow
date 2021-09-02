/************************************************************************
    > File Name: src/parser/tools.cpp
    > Author: deqi
    > Mail: deqi_tang@163.com 
    > Created Time: Tue 23 Feb 2021 09:38:15 PM CST
************************************************************************/

#include "askit/parser/tools.h"

#include <boost/filesystem.hpp>
#include "askit/parser/cif.h"
#include "askit/parser/xyz.h"
#include "askit/base/crystal.h"
#include "askit/utils.h"


namespace askit {
    
    namespace filesys = boost::filesystem;

    Crystal read_structure_file(std::string filepath) {
        // read structure file
        Crystal crystal;
        filesys::path in_path(filepath);
        //std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;

        if (in_path.has_extension() && in_path.extension().string() == ".cif") {
            askit::read_cif_file(&crystal, filepath);
        } else if (in_path.has_extension() && in_path.extension().string() == ".xyz") {
            askit::read_xyz_file(&crystal, filepath);
        }
        
        return crystal;
    }

} // namespace askit
