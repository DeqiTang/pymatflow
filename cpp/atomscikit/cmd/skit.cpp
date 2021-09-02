#include <boost/program_options.hpp>
//#include <experimental/filesystem>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <regex>
#include <cstdlib>
#include "askit/parser/cif.h"
#include "askit/parser/xyz.h"
#include "askit/base/crystal.h"
#include "askit/utils.h"

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
    std::cout << "*                       skit.x utils runnig                          *" << std::endl;
    std::cout << "**********************************************************************" << std::endl;
    
    if (0 == vm.count("command")) { // or by vm.empty()
        std::cout << "You didn't specify any subcommand!\n";
        std::exit(1);
    }
    std::string cmd = vm["command"].as<std::string>();


    if (cmd == "convert") {
        //
        std::cout << "----------------------------------------------------------------------" << std::endl;    
        std::cout << "sub commands -> convert                                               " << std::endl;
        std::cout << "Run info:" << std::endl;

        
        // convert command has the following options:
        po::options_description opt_convert("convert options");
        opt_convert.add_options()
            ("input, i", po::value<std::string>()->required(), "input structure file")
            ("output, o", po::value<std::string>()->required(), "output structure file");
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_convert).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_convert << std::endl;
            std::exit(1);
        }
        po::notify(vm);


        if (vm.count("input") && vm.count("output")) {
            askit::Crystal crystal;
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

            if (in_path.has_extension() && in_path.extension().string() == ".cif") {
                askit::read_cif_file(&crystal, input_file);
            } else if (in_path.has_extension() && in_path.extension().string() == ".xyz") {
                askit::read_xyz_file(&crystal, input_file);
            }
        
            // write structure file
        
            if (out_path.has_extension() && out_path.extension().string() == ".cif") {
                askit::write_cif_file(&crystal, output_file);
            } else if (out_path.has_extension() && out_path.extension().string() == ".xyz") {
                askit::write_xyz_file(&crystal, output_file);
            }
        }
        
        
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "sub command: convert finished!!!                                      " << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;
        
    } else if (cmd == "supercell") {
        //
        std::cout << "----------------------------------------------------------------------" << std::endl;    
        std::cout << "sub commands -> supercell                                             " << std::endl;
        std::cout << "Run info:" << std::endl;
        
        // supercell command has the following options:
        po::options_description opt_supercell("supercell options");
        opt_supercell.add_options()
            ("input, i", po::value<std::string>()->required(), "input structure file")
            ("output, o", po::value<std::string>()->required(), "output structure file")
            ("n", po::value<std::vector<int> >()->multitoken(), "n to define the supercell")
            ("precision", po::value<double>()->default_value(1.0e-8), "precision for redefine lattice");            
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin())        ;
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_supercell).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_supercell << std::endl;
            std::exit(1);
        }
        po::notify(vm);


        if (vm.count("input") && vm.count("output")) {
            askit::Crystal crystal;
        
            std::string input_file = vm["input"].as<std::string>();
            std::string output_file = vm["output"].as<std::string>();
            std::vector<int> n = vm["n"].as<std::vector<int>>();
            
            filesys::path in_path(input_file);
            filesys::path out_path(output_file);
        
            //std::cout << "input file: " << input_file << std::endl;
            //std::cout << "output file: " << output_file << std::endl; 
    
            //std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;
        
            // read structure file
            if (in_path.has_extension() && in_path.extension().string() == ".cif") {
                askit::read_cif_file(&crystal, input_file);
            } else if (in_path.has_extension() && in_path.extension().string() == ".xyz") {
                askit::read_xyz_file(&crystal, input_file);
            }
                        
            // build the supercell
            std::cout << "we are building the supercell from" << std::endl;
            std::cout << input_file << std::endl;
            std::cout << "to" << std::endl;
            std::cout << output_file << std::endl;
            std::cout << "by -> " << n[0] << "x" << n[1] << "x" << n[2] << std::endl;
            crystal.build_supercell(n);
                    
            // write structure file
        
            if (out_path.has_extension() && out_path.extension().string() == ".cif") {
                askit::write_cif_file(&crystal, output_file);
            } else if (out_path.has_extension() && out_path.extension().string() == ".xyz") {
                askit::write_xyz_file(&crystal, output_file);
            }
        }        
        // end supercell cmd
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "sub command: supercell finished!!!                                    " << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;
    } else if (cmd == "redefine") {
        //
        //
        std::cout << "----------------------------------------------------------------------" << std::endl;    
        std::cout << "sub commands -> redefine                                              " << std::endl;
        std::cout << "Run info:" << std::endl;
        
        // redefine command has the following options:
        po::options_description opt_redefine("redefine options");
        opt_redefine.add_options()
            ("input, i", po::value<std::string>(), "input structure file")
            ("output, o", po::value<std::string>(), "output structure file")
            ("a", po::value<std::vector<int> >()->multitoken(), "three integer to define the new a vector")
            ("b", po::value<std::vector<int> >()->multitoken(), "three integer to define the new b vector")
            ("c", po::value<std::vector<int> >()->multitoken(), "three integer to define the new c vector")
            ("precision", po::value<double>()->default_value(1.0e-8), "precision for redefine lattice");
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        //po::store(po::command_line_parser(opts).options(opt_redefine).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).run(), vm);
        po::store(po::command_line_parser(opts).options(opt_redefine).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_redefine << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        std::cout << "abc info" << std::endl;
        std::cout << "a: " << vm["a"].as<std::vector<int>>()[0] << " " << vm["a"].as<std::vector<int>>()[1] << " " << vm["a"].as<std::vector<int>>()[2] << std::endl;
        std::cout << "b: " << vm["b"].as<std::vector<int>>()[0] << " " << vm["b"].as<std::vector<int>>()[1] << " " << vm["b"].as<std::vector<int>>()[2] << std::endl;
        std::cout << "c: " << vm["c"].as<std::vector<int>>()[0] << " " << vm["c"].as<std::vector<int>>()[1] << " " << vm["c"].as<std::vector<int>>()[2] << std::endl;
        

        if (vm.count("input") && vm.count("output")) {
            askit::Crystal crystal;
        
        
            std::string input_file = vm["input"].as<std::string>();
            std::string output_file = vm["output"].as<std::string>();
            
            filesys::path in_path(input_file);
            filesys::path out_path(output_file);
        
        
            std::cout << "input: " << input_file << std::endl;
            std::cout << "output: " << output_file << std::endl;
        
    
        
            // read structure file
            std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;

            if (in_path.has_extension() && in_path.extension().string() == ".cif") {
                askit::read_cif_file(&crystal, input_file);
            } else if (in_path.has_extension() && in_path.extension().string() == ".xyz") {
                askit::read_xyz_file(&crystal, input_file);
            }
            
            // build the supercell
            std::vector<int> a = vm["a"].as<std::vector<int>>();
            std::vector<int> b = vm["b"].as<std::vector<int>>();
            std::vector<int> c = vm["c"].as<std::vector<int>>();

            
            redefine_lattice(&crystal, a, b, c, vm["precision"].as<double>());
            
        
            // write structure file
        
            if (out_path.has_extension() && out_path.extension().string() == ".cif") {
                askit::write_cif_file(&crystal, output_file);
            } else if (out_path.has_extension() && out_path.extension().string() == ".xyz") {
                askit::write_xyz_file(&crystal, output_file);
            }
        }        
        // end redefine lattice cmd
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "sub command: redefine finished!!!                                     " << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;        
    } else if (cmd == "tube") {
        //
        //
        std::cout << "----------------------------------------------------------------------" << std::endl;    
        std::cout << "sub commands -> tube                                                  " << std::endl;
        std::cout << "Run info:" << std::endl;
        
        // redefine command has the following options:
        po::options_description opt_tube("nanotube options");
        opt_tube.add_options()
            ("input, i", po::value<std::string>(), "input structure file")
            ("output, o", po::value<std::string>(), "output structure file")
            ("axis", po::value<std::string>(), "build nanotube along an axis parallel to axis specified, can be: a, b or c")
            ("plane", po::value<int>(), "on which plane to add vacuum layer. 1: ab, 2: ac, 3: bc");
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_tube).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_tube << std::endl;
            std::exit(1);
        }
        po::notify(vm);



        if (vm.count("input") && vm.count("output")) {
            askit::Crystal crystal;
        
        
            std::string input_file = vm["input"].as<std::string>();
            std::string output_file = vm["output"].as<std::string>();
            
            filesys::path in_path(input_file);
            filesys::path out_path(output_file);
        
        
            std::cout << "input: " << input_file << std::endl;
            std::cout << "output: " << output_file << std::endl;
        
    
        
            // read structure file
            std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;

            if (in_path.has_extension() && in_path.extension().string() == ".cif") {
                askit::read_cif_file(&crystal, input_file);
            } else if (in_path.has_extension() && in_path.extension().string() == ".xyz") {
                askit::read_xyz_file(&crystal, input_file);
            }
            
    
            // build the nanotube

            if (vm["plane"].as<int>() == 1) {
                crystal = build_nanotube_ab(&crystal, vm["axis"].as<std::string>());
            } else if (vm["plane"].as<int>() == 2) {
                crystal = build_nanotube_ac(&crystal, vm["axis"].as<std::string>());
            } else if (vm["plane"].as<int>() == 3) {
                crystal = build_nanotube_bc(&crystal, vm["axis"].as<std::string>());
            }
            
        
            // write structure file
        
            if (out_path.has_extension() && out_path.extension().string() == ".cif") {
                askit::write_cif_file(&crystal, output_file);
            } else if (out_path.has_extension() && out_path.extension().string() == ".xyz") {
                askit::write_xyz_file(&crystal, output_file);
            }
        }        
        // end supercell cmd
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "sub command: tube finished!!!                                         " << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;           
        
    } else if (cmd == "merge") {
        //
        //
        std::cout << "----------------------------------------------------------------------" << std::endl;    
        std::cout << "sub commands -> merge                                                 " << std::endl;
        std::cout << "Run info:" << std::endl;
        
        // redefine command has the following options:
        po::options_description opt_merge("merge layers options");
        opt_merge.add_options()
            ("input, i", po::value<std::vector<std::string>>()->multitoken()->required(), "input structure file")
            ("output, o", po::value<std::string>()->required(), "output structure file")
            ("usecell", po::value<int>()->default_value(0), "use cell of structure 1 or 2, otherwise 0 (avearge) by default")
            ("thick", po::value<double>()->default_value(10), "thickness of the vacuum layer")
            ("distance", po::value<double>()->default_value(3.4), "distance between the layer, in unit of Angstromg, default is 3.4");
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_merge).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_merge << std::endl;
            std::exit(1);
        }
        po::notify(vm);



        if (vm.count("input") && vm.count("output")) {
            askit::Crystal crystal1;
            askit::Crystal crystal2;
            askit::Crystal crystal;
        
        
            std::vector<std::string> input_files = vm["input"].as<std::vector<std::string>>();
            std::string output_file = vm["output"].as<std::string>();
            
            std::vector<filesys::path> in_paths;
            in_paths.push_back(filesys::path{input_files[0]});
            in_paths.push_back(filesys::path{input_files[1]});
            filesys::path out_path(output_file);
        
        
            std::cout << "input: " << input_files[0] << " " << input_files[1] << std::endl;
            std::cout << "output: " << output_file << std::endl;
        
    
        
            // read structure file
            std::cout << "in_path(extension)->" << in_paths[0].extension().string() << std::endl;
            std::cout << "in_path(extension)->" << in_paths[1].extension().string() << std::endl;

            if (in_paths[0].has_extension() && in_paths[0].extension().string() == ".cif") {
                askit::read_cif_file(&crystal1, input_files[0]);
            } else if (in_paths[0].has_extension() && in_paths[0].extension().string() == ".xyz") {
                askit::read_xyz_file(&crystal1, input_files[0]);
            }
            
            if (in_paths[1].has_extension() && in_paths[1].extension().string() == ".cif") {
                askit::read_cif_file(&crystal2, input_files[1]);
            } else if (in_paths[1].has_extension() && in_paths[1].extension().string() == ".xyz") {
                askit::read_xyz_file(&crystal2, input_files[1]);
            }
            
    
            // merge the layers
            std::cout << "beginning to merge layers" << std::endl;
            crystal = merge_layers(&crystal1, &crystal2, vm["usecell"].as<int>(), vm["distance"].as<double>(), vm["thick"].as<double>());
        
            // write structure file
        
            if (out_path.has_extension() && out_path.extension().string() == ".cif") {
                askit::write_cif_file(&crystal, output_file);
            } else if (out_path.has_extension() && out_path.extension().string() == ".xyz") {
                askit::write_xyz_file(&crystal, output_file);
            }
        }        
        // end merge layers cmd
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "sub command: merge finished!!!                                        " << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;               
    } else if (cmd == "cleave") {
        
        //
        //
        std::cout << "----------------------------------------------------------------------" << std::endl;    
        std::cout << "sub commands -> cleave                                                " << std::endl;
        std::cout << "Run info:" << std::endl;
        
        // redefine command has the following options:
        po::options_description opt_cleave("cleave options");
        opt_cleave.add_options()
            ("input, i", po::value<std::string>(), "input structure file")
            ("output, o", po::value<std::string>(), "output structure file")
            ("direction", po::value<std::vector<int>>()->multitoken(), "directino of the surface plane to cleave")
            ("thick", po::value<double>()->default_value(10), "thickness of the vacuum layer")
            ("precision", po::value<double>()->default_value(1.0e-8), "a value that is large than 0 and infinitely close to 0 used to judge whether one atom is in another periodic of the redefined cell used in cleave surface");
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_cleave).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_cleave << std::endl;
            std::exit(1);
        }
        po::notify(vm);



        if (vm.count("input") && vm.count("output")) {
            askit::Crystal crystal;
        
        
            std::string input_file = vm["input"].as<std::string>();
            std::string output_file = vm["output"].as<std::string>();
            
            filesys::path in_path(input_file);
            filesys::path out_path(output_file);
        
        
            std::cout << "input: " << input_file << std::endl;
            std::cout << "output: " << output_file << std::endl;
        
    
        
            // read structure file
            std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;

            if (in_path.has_extension() && in_path.extension().string() == ".cif") {
                askit::read_cif_file(&crystal, input_file);
            } else if (in_path.has_extension() && in_path.extension().string() == ".xyz") {
                askit::read_xyz_file(&crystal, input_file);
            }
            
            std::cout << "finish reading the structure file, begin to cleave the surface" << std::endl;
    
            // cleave surface
            cleave_surface(&crystal, vm["direction"].as<std::vector<int>>(), vm["thick"].as<double>(), vm["precision"].as<double>());
            
            std::cout << "finish cleaving the surface" << std::endl;
            // write structure file
        
            if (out_path.has_extension() && out_path.extension().string() == ".cif") {
                askit::write_cif_file(&crystal, output_file);
            } else if (out_path.has_extension() && out_path.extension().string() == ".xyz") {
                askit::write_xyz_file(&crystal, output_file);
            }
        }        
        // end cleave surface cmd
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "sub command: cleave finished!!!                                       " << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;            
                
    } else if (cmd == "vacuum") {

        //
        std::cout << "----------------------------------------------------------------------" << std::endl;    
        std::cout << "sub commands -> vacuum                                                " << std::endl;
        std::cout << "Run info:" << std::endl;
        
        // redefine command has the following options:
        po::options_description opt_vacuum("vacuum options");
        opt_vacuum.add_options()
            ("input, i", po::value<std::string>(), "input structure file")
            ("output, o", po::value<std::string>(), "output structure file")
            ("plane", po::value<int>(), "plane to add vacuum layer, 1: ab, 2: ac, 3: bc")
            ("thick", po::value<double>()->default_value(10), "thickness of the vacuum layer");
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_vacuum).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);

        if (vm.count("help")) {
            std::cout << opt_vacuum << std::endl;
            std::exit(1);
        }
        po::notify(vm);



        if (vm.count("input") && vm.count("output")) {
            askit::Crystal crystal;
        
        
            std::string input_file = vm["input"].as<std::string>();
            std::string output_file = vm["output"].as<std::string>();
            
            filesys::path in_path(input_file);
            filesys::path out_path(output_file);
        
        
            std::cout << "input: " << input_file << std::endl;
            std::cout << "output: " << output_file << std::endl;
        
    
        
            // read structure file
            std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;

            if (in_path.has_extension() && in_path.extension().string() == ".cif") {
                askit::read_cif_file(&crystal, input_file);
            } else if (in_path.has_extension() && in_path.extension().string() == ".xyz") {
                askit::read_xyz_file(&crystal, input_file);
            }
            
            std::cout << "finish reading the structure file, begin to add vacuum layer" << std::endl;
    
            // add vacuum layer
            vacuum_layer(&crystal, vm["plane"].as<int>(), vm["thick"].as<double>());
            
            std::cout << "finish add vacuum layer" << std::endl;
            // write structure file
        
            if (out_path.has_extension() && out_path.extension().string() == ".cif") {
                askit::write_cif_file(&crystal, output_file);
            } else if (out_path.has_extension() && out_path.extension().string() == ".xyz") {
                askit::write_xyz_file(&crystal, output_file);
            }
        }        
        // end supercell cmd
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "sub command: vacuum  finished!!!                                      " << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;   
        
    } else {
        std::cout << "The specified subcommand is not defined!\n";
    }

    //
    return 0;
}
