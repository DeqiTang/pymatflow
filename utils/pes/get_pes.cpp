/*
get_pes.sh running too slowly. this get_pes.cpp can help with it.
*/

#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <regex>
#include <boost/algorithm/string.hpp>
#include <cstdlib>

namespace fs = std::filesystem;


double get_final_energy_from_outcar(std::ifstream& outcar) {
    //
    std::regex energy_regex(".+entropy=.+");
    std::string line;
    //std::vector<std::string> energy_strings;
    std::string energy_string;
    double final_energy;

    while (std::getline(outcar, line)) {
        if (std::regex_match(line, energy_regex)) {
            //energy_strings.push_back(line);
            energy_string = line;
        }
    }

    // using boost:split to split the string
    //std::vector<std::string> energy_str_split;
    //boost::split(energy_str_split, energy_string, boost::is_any_of(" "));

    
    //also we can use std::regex to implement it
    std::regex whitespace("\\s+");
    std::vector<std::string> energy_str_split(std::sregex_token_iterator(energy_string.begin(), energy_string.end(), whitespace, -1), 
        std::sregex_token_iterator());

    //return double(energy_str_split[3]);
    std::cout << energy_str_split[4] << std::endl;
    return (double)std::atof(energy_str_split[4].c_str());
};


int main() {

    std::vector<std::string> xy_dirs;
    // regex: _numeric_numeric_
    std::regex xy_dir_regex("_(-[0-9]+(.[0-9]+)?|[0-9]+(.[0-9]+)?)_(-[0-9]+(.[0-9]+)?|[0-9]+(.[0-9]+)?)_");
    
    for (auto& p : fs::directory_iterator("./")) {
        fs::path path(p.path());
        //std::cout << path.filename().string() << std::endl;
        if (p.is_directory() && std::regex_match(path.filename().string(), xy_dir_regex)) {
            xy_dirs.push_back(path.filename().string());
            std::cout << path.filename().string() << std::endl;
        }
    }

    std::ifstream opt_out;
    std::ofstream pes_data_file;
    fs::create_directory("post-processing");
    pes_data_file.open("post-processing/pes.data", std::ios::out);
    pes_data_file << "#format: x y energy\n";

    for (auto dir : xy_dirs) {
        opt_out.open(dir+"/OUTCAR");
        /*&
        while (std::getline(opt_out, line)) {
            if (std::regex_match(line, energy_regex)) {
                energy_strings.push_back(line);
            }
        }
        */
        //std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1) << get_final_energy_from_outcar(opt_out) << std::endl;
        std::regex xy_dir("_");
        std::vector<std::string> xy(std::sregex_token_iterator(dir.begin(), dir.end(), xy_dir, -1), 
            std::sregex_token_iterator());
        pes_data_file<< std::setprecision(std::numeric_limits<double>::digits10 + 1)  << xy[1] << " " << xy[2] << " " << get_final_energy_from_outcar(opt_out) << "\n";
    }

    return 0;
}