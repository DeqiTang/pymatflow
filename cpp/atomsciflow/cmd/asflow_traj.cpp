#include <boost/program_options.hpp>
//#include <experimental/filesystem>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <regex>
#include <cstdlib>
#include <armadillo>
#include "atomsciflow/parser/cif.h"
#include "atomsciflow/parser/xyz.h"
#include "atomsciflow/base/crystal.h"
#include "atomsciflow/parser/xyztraj.h"
#include "atomsciflow/utils.h"

// needs: libboost-dev, libboost-program-options-dev

namespace po = boost::program_options;


//namespace filesys = std::experimental::filesystem;  // --std=c++17 -lstdc++fs
namespace filesys = boost::filesystem;     // --std=c++11 -lboost_filesystem -lboost_system



// used to allow negative number parsed to boost cmd option
std::vector<po::option> ignore_numbers(std::vector<std::string>& args)
{
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
    std::cout << "*                       skit-traj.x utils runnig                     *" << std::endl;
    std::cout << "**********************************************************************" << std::endl;
    
    if (0 == vm.count("command")) { // or by vm.empty()
        std::cout << "You didn't specify any subcommand!\n";
        std::exit(1);
    }
    std::string cmd = vm["command"].as<std::string>();


    if (cmd == "tidy") {
        //
        std::cout << "----------------------------------------------------------------------" << std::endl;    
        std::cout << "sub commands -> tidy                                               " << std::endl;
        std::cout << "Run info:" << std::endl;

        
        // tidy command has the following options:
        po::options_description opt_tidy("tidy options");
        opt_tidy.add_options()
            ("input, i", po::value<std::string>()->required(), "input xyz trajectory file")
            ("output, o", po::value<std::string>()->required(), "output xyz trajectory file")
            ("cell", po::value<std::vector<double>>()->multitoken()->required(), "cell parameter");
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_tidy).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_tidy << std::endl;
            std::exit(1);
        }
        po::notify(vm);


        if (vm.count("input") && vm.count("output")) {
            std::vector<atomsciflow::Crystal> trajectory;
            //crystal.read_xyz_file(vm["input"].as<std::string>());
            //crystal.write_cif_file(vm["output"].as<std::string>());
        
        
            std::string input_file = vm["input"].as<std::string>();
            std::string output_file = vm["output"].as<std::string>();
        
            filesys::path in_path(input_file);
            filesys::path out_path(output_file);
        
        
            std::cout << "input: " << input_file << std::endl;
            std::cout << "output: " << output_file << std::endl;
        
            // read structure file
            std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;

            
            atomsciflow::read_xyztraj_file(trajectory, input_file, 0);

            std::cout << "num of images: " << trajectory.size() << std::endl;

            std::vector<std::vector<double>> cell;
            std::vector<double> arg_cell = vm["cell"].as<std::vector<double>>();
            std::vector<double> vec;
            vec.push_back(arg_cell[0]);
            vec.push_back(arg_cell[1]);
            vec.push_back(arg_cell[2]);
            cell.push_back(vec);
            vec.clear();
            vec.push_back(arg_cell[3]);
            vec.push_back(arg_cell[4]);
            vec.push_back(arg_cell[5]);
            cell.push_back(vec);
            vec.clear();
            vec.push_back(arg_cell[6]);
            vec.push_back(arg_cell[7]);
            vec.push_back(arg_cell[8]);
            cell.push_back(vec);
            vec.clear();
        
            std::cout << "cell processed" << std::endl;

            for (auto& image : trajectory) {
                image.cell = cell;
                set_frac_within_zero_and_one(&image);
            }
            std::cout << "tidied" << std::endl;
            // write structure file
            atomsciflow::write_xyztraj_file(trajectory, output_file, 0);
        }
        
        
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "sub command: tidy finished!!!                                      " << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;
        
    } else {
        std::cout << "The specified subcommand is not defined!\n";
    }
    
    //
    return 0;
}
