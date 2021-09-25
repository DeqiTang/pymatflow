#include <boost/program_options.hpp>
//#include <experimental/filesystem>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <regex>
#include <cstdlib>
#include "atomsciflow/parser/cif.h"
#include "atomsciflow/parser/xyz.h"
#include "atomsciflow/parser/tools.h"
#include "atomsciflow/base/crystal.h"
#include "atomsciflow/utils.h"
#include "atomsciflow/abinit/abinit.h"
#include "atomsciflow/cp2k/cp2k.h"
#include "atomsciflow/qe/pw.h"

// needs: libboost-dev, libboost-program-options-dev

namespace po = boost::program_options;


//namespace filesys = std::experimental::filesystem;  // --std=c++17 -lstdc++fs
namespace filesys = boost::filesystem;     // --std=c++11 -lboost_filesystem -lboost_system

int log_sub_cmd_start(std::string cmd) {
    std::cout << "----------------------------------------------------------------------" << std::endl;    
    std::cout << "sub commands -> " << cmd << std::endl;
    std::cout << "Run info:" << std::endl;
    return 0;
}

int log_sub_cmd_end(std::string cmd) {
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "sub command: " << cmd << " finished!!!                                      " << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
    return 0;
}

// used to allow negative number parsed to boost cmd option
std::vector<po::option> ignore_numbers(std::vector<std::string>& args) {
    // this function can help to alow negative number args but it probhibits positional args
    // however we do not need positional args. so it is ok.
    std::vector<po::option> result;
    int pos = 0;
    while(!args.empty()) {
        const auto& arg = args[0];
        double num;
        if(boost::conversion::try_lexical_convert(arg, num)) {
            result.push_back(po::option());
            po::option& opt = result.back();

            opt.position_key = pos++;
            opt.value.push_back(arg);
            opt.original_tokens.push_back(arg);

            args.erase(args.begin());
        } else {
            break;
        }
    }

    return result;
}



int main(int argc, char const* argv[]) {
    //

    po::options_description global("Global options");
    global.add_options()
        ("command", po::value<std::string>(), "command to execute")
        ("subargs", po::value<std::vector<std::string> >(), "Arguments for command")
        ("help, h", "print out help information");
        
    po::positional_options_description pos;
    pos.add("command", 1).add("subargs", -1);
    
    po::variables_map vm;
    
    po::parsed_options parsed = po::command_line_parser(argc, argv).
        options(global).
        positional(pos).
        allow_unregistered().
        run();
        
    po::store(parsed, vm);

    std::cout << "**********************************************************************" << std::endl;
    std::cout << "*                       askit.x utils runnig                         *" << std::endl;
    std::cout << "**********************************************************************" << std::endl;
    
    if (0 == vm.count("command")) { // or by vm.empty()
        std::cout << "You didn't specify any subcommand!\n";
        std::exit(1);
    }
    std::string cmd = vm["command"].as<std::string>();
    

    if (cmd == "abinit") {
        //
        log_sub_cmd_start(cmd);
        
        // convert command has the following options:
        po::options_description opt_abinit("abinit options");
        opt_abinit.add_options()
            ("input, i", po::value<std::string>()->required(), "input structure file");
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_abinit).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_abinit << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        if (vm.count("input")) {
            atomsciflow::Crystal crystal;
            //crystal.read_xyz_file(vm["input"].as<std::string>());
            //crystal.write_cif_file(vm["output"].as<std::string>());
        
        
            std::string input_file = vm["input"].as<std::string>();
            //std::string output_file = vm["output"].as<std::string>();
        
            filesys::path in_path(input_file);
            //filesys::path out_path(output_file);
        
        
            std::cout << "input: " << input_file << std::endl;
            //std::cout << "output: " << output_file << std::endl;
        
    
        
            // read structure file
            //std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;
            crystal = atomsciflow::read_structure_file(input_file);
            
            atomsciflow::Abinit calculator;
            std::cout << calculator.to_string() << std::endl;
        }
        
        log_sub_cmd_end(cmd); 
        
    } else if (cmd == "cp2k") {
        //
        log_sub_cmd_start(cmd);
        
        // convert command has the following options:
        po::options_description opt_cp2k("cp2k options");
        opt_cp2k.add_options()
            ("input, i", po::value<std::string>()->required(), "input structure file")
            ("directory, d", po::value<std::string>()->default_value("askit-calc-running"));
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_cp2k).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_cp2k << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        if (vm.count("input")) {
            atomsciflow::Crystal crystal;
            //crystal.read_xyz_file(vm["input"].as<std::string>());
            //crystal.write_cif_file(vm["output"].as<std::string>());
        
        
            std::string input_file = vm["input"].as<std::string>();
            //std::string output_file = vm["output"].as<std::string>();
        
            filesys::path in_path(input_file);
            //filesys::path out_path(output_file);
        
        
            std::cout << "input: " << input_file << std::endl;
            //std::cout << "output: " << output_file << std::endl;
        
            // read structure file
            // std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;
            crystal = atomsciflow::read_structure_file(input_file); 

            atomsciflow::Cp2k calculator;
            calculator.set_subsys(crystal);
            std::cout << calculator.to_string() << std::endl;

        }
       
        log_sub_cmd_end(cmd);
        
    } else if ("qe" == cmd) {
        //
        log_sub_cmd_start(cmd);
        
        // convert command has the following options:
        po::options_description opt_qe("qe options");
        opt_qe.add_options()
            ("input, i", po::value<std::string>()->required(), "input structure file")
            ("directory, d", po::value<std::string>()->default_value("askit-calc-running"));
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_qe).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_qe << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        if (vm.count("input")) {
            atomsciflow::Crystal crystal;
            //crystal.read_xyz_file(vm["input"].as<std::string>());
            //crystal.write_cif_file(vm["output"].as<std::string>());
        
        
            std::string input_file = vm["input"].as<std::string>();
            //std::string output_file = vm["output"].as<std::string>();
        
            filesys::path in_path(input_file);
            //filesys::path out_path(output_file);
        
        
            std::cout << "input: " << input_file << std::endl;
            //std::cout << "output: " << output_file << std::endl;
        
    
            // read structure file
            // std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;
            crystal = atomsciflow::read_structure_file(input_file);

            atomsciflow::QePw calculator;
            //calculator.set_subsys(crystal);
            std::cout << calculator.to_string() << std::endl;

        }
        
        log_sub_cmd_end(cmd);    
    } else {
        std::cout << "The specified subcommand is not defined!\n";
    }
        
    //
    return 0;
}
