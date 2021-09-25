#include <boost/program_options.hpp>
//#include <experimental/filesystem>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <algorithm>
#include <regex>
#include <cstdlib>
#include <cmath>
#include <armadillo>

#include "atomsciflow/parser/cif.h"
#include "atomsciflow/parser/xyz.h"
#include "atomsciflow/base/crystal.h"
#include "atomsciflow/utils.h"

#include "atomsciflow/cube_handle/cube_handle.h"
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
    std::cout << "*                       skit-cube.x utils runnig                     *" << std::endl;
    std::cout << "**********************************************************************" << std::endl;
    
    if (0 == vm.count("command")) { // or by vm.empty()
        std::cout << "You didn't specify any subcommand!\n";
        std::exit(1);
    }
    std::string cmd = vm["command"].as<std::string>();


    if (cmd == "along") {
        //
        std::cout << "----------------------------------------------------------------------" << std::endl;    
        std::cout << "sub commands -> along                                                 " << std::endl;
        std::cout << "Run info:" << std::endl;

        // along command has the following options:
        po::options_description opt_along("along options");
        opt_along.add_options()
            ("input, i", po::value<std::string>(), "input cube file")
            ("output, o", po::value<std::string>(), "output structure file")
            ("abscissa", po::value<std::vector<std::string> >()->multitoken(), "choose the direction to do the dimension reduction");
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_along).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);
        std::cout << "parse arguments again" << "\n";

        if (vm.count("help")) {
            std::cout << opt_along << std::endl;
            std::exit(1);
        }
        po::notify(vm);



        std::cout << "Preparing to deal with input and output\n";

        if (vm.count("input") && vm.count("output")) {

        
            std::string input_file = vm["input"].as<std::string>();
            std::string output_file = vm["output"].as<std::string>();

            filesys::path in_path(input_file);
            filesys::path out_path(output_file);

            std::vector<std::string> abscissa = vm["abscissa"].as<std::vector<std::string>>();
        
            std::cout << "input: " << input_file << std::endl;
            std::cout << "output: " << output_file << std::endl;
        
    
        
            // read structure file
            std::cout << "in_path(extension)->" << in_path.extension().string() << std::endl;

            
            atomsciflow::CubeElectronDensity cube;
            cube.read_cube(input_file);
        
            double bohr_to_angstrom = 0.529177249;
            // data dimension reduction
            // the unit of value is actually not physical now!
            // cell_volume are in unit of Angstrom^3
            arma::mat latcell(3, 3);
            latcell.row(0) = arma::conv_to<arma::rowvec>::from(cube.crystal.cell[0]);
            latcell.row(1) = arma::conv_to<arma::rowvec>::from(cube.crystal.cell[1]);
            latcell.row(2) = arma::conv_to<arma::rowvec>::from(cube.crystal.cell[2]);
            
            std::cout << "latcell:\n";
            std::cout << latcell << std::endl;      

            // std::cout << "test value:\n";
            // std::cout << arma::dot(arma::cross(latcell.row(0), latcell.row(1)), latcell.row(2)) << std::endl;

            double cell_volume = arma::dot(arma::cross(latcell.row(0), latcell.row(1)), latcell.row(2));

            double cell_volume_per_unit = cell_volume / (cube.ngridx * cube.ngridy * cube.ngridz);

            //std::cout << "cell volume per unit: " << cell_volume_per_unit << std::endl;

            int ngridx = cube.ngridx;
            int ngridy = cube.ngridy;
            int ngridz = cube.ngridz;
            // value in cube file are \rho(r)_of_electrons in unit of e/Bohr^3
            // namely number of electrons each Borh^3
            // so we have to convert it to e/Angstrom^3, through divide it by borh_to_angstrom**3
        
            double total_electrons = arma::accu(cube.data)  * cell_volume_per_unit / pow(bohr_to_angstrom, 3);

            std::cout << "======================================================\n";
            std::cout << "           Information collected\n";
            std::cout << "------------------------------------------------------\n";
            std::cout << "cell volume: " << cell_volume << " (A^3)\n";
            std::cout << "total electrons: " << total_electrons << "\n";
            
            // cube.data is in unit of e/Bohr^3
            // we will build data_red_? to be in unit of e/Anstrom, namely number of electrons per Angstrom
            // to do this we have to time the volume density with bohr_to_angstrom^-3
            std::vector<double> data_red_a;
            std::vector<double> data_red_b;
            std::vector<double> data_red_c;
            //arma::vec data_red_a;
            //arma::vec data_red_b;
            //arma::vec data_red_c;

            double a = arma::norm(latcell.row(0), 2);
            double b = arma::norm(latcell.row(1), 2);
            double c = arma::norm(latcell.row(2), 2);
            
            if (std::find(abscissa.begin(), abscissa.end(), "c") != abscissa.end()) {
                double len_ci = c / ngridz;
                auto size = arma::size(cube.data);
                double tmp;
                double rho_line;
                for (int ci = 0; ci < size[2]; ci++) {
                    tmp  = 0;
                    for (int bi = 0; bi < size[1]; bi++) {
                        for (int ai = 0; ai < size[0]; ai++) {
                            tmp += cube.data.at(ai, bi, ci);
                        }
                    }
                    rho_line = tmp * cell_volume_per_unit / pow(bohr_to_angstrom, 3) / len_ci;
                    data_red_c.push_back(rho_line);
                }
            }

            if (std::find(abscissa.begin(), abscissa.end(), "b") != abscissa.end())  {
                double len_bi = b / ngridy;
                auto size = arma::size(cube.data);
                double tmp;
                double rho_line;
                for (int bi  = 0; bi < size[1]; bi++) {
                    tmp = 0;
                    for (int ai = 0; ai < size[0]; ai++) {
                        for (int ci = 0; ci < size[2]; ci++) {
                            tmp += cube.data.at(ai, bi, ci);
                        }
                    }
                    rho_line = tmp * cell_volume_per_unit / pow(bohr_to_angstrom, 3) / len_bi;
                    data_red_b.push_back(rho_line);
                }
            }

            if (std::find(abscissa.begin(), abscissa.end(), "a") != abscissa.end())  {
                double len_ai = a / ngridx;
                auto size = arma::size(cube.data);
                double tmp;
                double rho_line;
                for (int ai = 0; ai < size[0]; ai++) {
                    tmp = 0;
                    for (int ci = 0; ci < size[2]; ci++) {
                        for (int bi = 0; bi < size[1]; bi++) {
                            tmp += cube.data.at(ai, bi, ci);
                        }
                    }
                    rho_line = tmp * cell_volume_per_unit / pow(bohr_to_angstrom, 3) / len_ai;
                    data_red_a.push_back(rho_line);
                }
            }

            std::cout << "Preparing to output the data\n";
            // output the data and make the plot
            if (std::find(abscissa.begin(), abscissa.end(), "c") != abscissa.end()) {
                std::ofstream fout;
                fout.open(output_file+".1d.c.data");
                fout << "#c(angstrom) rho(e) (number of electron per Angstrom)\n";
                arma::vec c_coord = arma::linspace(0, c, data_red_c.size());
                for (int i = 0; i < data_red_c.size(); i++) {
                    fout << c_coord.at(i) << " " << data_red_c.at(i) << "\n";
                }
                fout.close();
            }
           
            if (std::find(abscissa.begin(), abscissa.end(), "b") != abscissa.end()) {
                std::ofstream fout;
                fout.open(output_file+".1d.b.data");
                fout << "#b(angstrom) rho(e) (number of electron per Angstrom)\n";
                arma::vec b_coord = arma::linspace(0, b, data_red_b.size());
                for (int i = 0; i < data_red_b.size(); i++) {
                    fout << b_coord.at(i) << " " << data_red_b.at(i) << "\n";
                }
                fout.close();
            }

            if (std::find(abscissa.begin(), abscissa.end(), "a") != abscissa.end())  {
                std::ofstream fout;
                fout.open(output_file+".1d.a.data");
                fout << "#a(angstrom) rho(e) (number of electron per Angstrom)\n";
                arma::vec a_coord = arma::linspace(0, a, data_red_a.size());
                for (int i = 0; i < data_red_a.size(); i++) {
                    fout << a_coord.at(i) << " " << data_red_a.at(i) << "\n";
                }
                fout.close();
            }       

            //std::ofstream fxx;
            //fxx.open("xxx.data");
            //fxx << cube.data << std::endl;
            //fxx.close();

        }
        
        
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "sub command: along finished!!!                                      " << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;
        
    } else if (cmd == "diff") {
        //
        std::cout << "----------------------------------------------------------------------" << std::endl;    
        std::cout << "sub commands -> diff                                                 " << std::endl;
        std::cout << "Run info:" << std::endl;

        // along command has the following options:
        po::options_description opt_diff("diff options");
        opt_diff.add_options()
            ("input, i", po::value<std::vector<std::string> >()->multitoken(), "input three cube file: TOTAL PART1 PART2")
            ("output, o", po::value<std::string>(), "output structure file")
            ("abscissa", po::value<std::vector<std::string> >()->multitoken(), "choose the direction to do the dimension reduction");
        //opts.add_options()
        //    ("input, i", po::value<std::string>(&input), "input structure file")
        //    ("output, o", po::value<std::string>(&output), "output structure file");
        // collect all the unrecognized options from the first pass. this will include the 
        // (positional) command name so we need to erase that
        std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
        opts.erase(opts.begin());
        //parse again...
        po::store(po::command_line_parser(opts).options(opt_diff).style(po::command_line_style::unix_style | po::command_line_style::allow_long_disguise).extra_style_parser(&ignore_numbers).run(), vm);
        std::cout << "parse arguments again" << "\n";

        if (vm.count("help")) {
            std::cout << opt_diff << std::endl;
            std::exit(1);
        }
        po::notify(vm);



        std::cout << "Preparing to deal with input and output\n";

        if (vm.count("input") && vm.count("output")) {

        
            // std::string input_file = vm["input"].as<std::string>();
            std::vector<std::string> input_files = vm["input"].as<std::vector<std::string>>();
            std::string output_file = vm["output"].as<std::string>();

            // filesys::path in_path(input_file);
            filesys::path out_path(output_file);

            std::vector<std::string> abscissa = vm["abscissa"].as<std::vector<std::string>>();
        
            std::cout << "input: \n" 
                << input_files[0] << "\n" 
                << input_files[1] << "\n" 
                << input_files[2] << "\n" << std::endl;

            std::cout << "output: " << output_file << std::endl;
        

            atomsciflow::cube_diff_1d(input_files, output_file, abscissa);

        }
        
        
        std::cout << "----------------------------------------------------------------------" << std::endl;
        std::cout << "sub command: diff finished!!!                                      " << std::endl;
        std::cout << "----------------------------------------------------------------------" << std::endl;
        
    } else {
        std::cout << "The specified subcommand is not defined!\n";
    } 
    
    //
    return 0;
}
