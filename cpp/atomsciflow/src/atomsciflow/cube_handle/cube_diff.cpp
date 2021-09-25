#include "atomsciflow/cube_handle/cube_handle.h"

namespace atomsciflow {

int cube_diff_1d(std::vector<std::string> input_files, std::string output_file, std::vector<std::string> abscissa) {
    
    std::vector<atomsciflow::CubeElectronDensity> cubes;

    for (auto i = 0; i < 3; i++) {
        atomsciflow::CubeElectronDensity cube;
        cubes.push_back(cube);
    }
    cubes[0].read_cube(input_files[0]);
    cubes[1].read_cube(input_files[1]);
    cubes[2].read_cube(input_files[2]);

    double bohr_to_angstrom = 0.529177249;
    // data dimension reduction
    // the unit of value is actually not physical now!
    // cell_volume are in unit of Angstrom^3
    arma::mat latcell(3, 3);
    latcell.row(0) = arma::conv_to<arma::rowvec>::from(cubes[0].crystal.cell[0]);
    latcell.row(1) = arma::conv_to<arma::rowvec>::from(cubes[0].crystal.cell[1]);
    latcell.row(2) = arma::conv_to<arma::rowvec>::from(cubes[0].crystal.cell[2]);
    
    std::cout << "latcell:\n";
    std::cout << latcell << std::endl;   

    // std::cout << "test value:\n";
    // std::cout << arma::dot(arma::cross(latcell.row(0), latcell.row(1)), latcell.row(2)) << std::endl;

    double cell_volume = arma::dot(arma::cross(latcell.row(0), latcell.row(1)), latcell.row(2));

    double cell_volume_per_unit = cell_volume / (cubes[0].ngridx * cubes[0].ngridy * cubes[0].ngridz);

    //std::cout << "cell volume per unit: " << cell_volume_per_unit << std::endl;

    int ngridx = cubes[0].ngridx;
    int ngridy = cubes[0].ngridy;
    int ngridz = cubes[0].ngridz;
    // value in cube file are \rho(r)_of_electrons in unit of e/Bohr^3
    // namely number of electrons each Borh^3
    // so we have to convert it to e/Angstrom^3, through divide it by borh_to_angstrom**3

    double total_electrons = arma::accu(cubes[0].data)  * cell_volume_per_unit / pow(bohr_to_angstrom, 3);

    std::cout << "======================================================\n";
    std::cout << "           Information collected\n";
    std::cout << "------------------------------------------------------\n";
    std::cout << "cell volume: " << cell_volume << " (A^3)\n";
    std::cout << "total electrons: " << total_electrons << "\n";
    
    // cube.data is in unit of e/Bohr^3
    // we will build data_red_? to be in unit of e/Anstrom, namely number of electrons per Angstrom
    // to do this we have to time the volume density with bohr_to_angstrom^-3

    arma::cube diff_data = cubes[0].data - cubes[1].data - cubes[2].data;

    std::vector<double> diff_data_red_a;
    std::vector<double> diff_data_red_b;
    std::vector<double> diff_data_red_c;
    //arma::vec diff_data_red_a;
    //arma::vec diff_data_red_b;
    //arma::vec diff_data_red_c;

    double a = arma::norm(latcell.row(0), 2);
    double b = arma::norm(latcell.row(1), 2);
    double c = arma::norm(latcell.row(2), 2);
    
    auto size = arma::size(cubes[0].data);

    if (std::find(abscissa.begin(), abscissa.end(), "c") != abscissa.end()) {
        double len_ci = c / ngridz;
        //auto size = arma::size(cubes[0].data);
        double tmp;
        double rho_line;
        for (int ci = 0; ci < size[2]; ci++) {
            tmp  = 0;
            for (int bi = 0; bi < size[1]; bi++) {
                for (int ai = 0; ai < size[0]; ai++) {
                    tmp += diff_data.at(ai, bi, ci);
                }
            }
            rho_line = tmp * cell_volume_per_unit / pow(bohr_to_angstrom, 3) / len_ci;
            diff_data_red_c.push_back(rho_line);
        }
    }

    if (std::find(abscissa.begin(), abscissa.end(), "b") != abscissa.end())  {
        double len_bi = b / ngridy;
        // auto size = arma::size(cube.data);
        double tmp;
        double rho_line;
        for (int bi  = 0; bi < size[1]; bi++) {
            tmp = 0;
            for (int ai = 0; ai < size[0]; ai++) {
                for (int ci = 0; ci < size[2]; ci++) {
                    tmp += diff_data.at(ai, bi, ci);
                }
            }
            rho_line = tmp * cell_volume_per_unit / pow(bohr_to_angstrom, 3) / len_bi;
            diff_data_red_b.push_back(rho_line);
        }
    }

    if (std::find(abscissa.begin(), abscissa.end(), "a") != abscissa.end())  {
        double len_ai = a / ngridx;
        // auto size = arma::size(cube.data);
        double tmp;
        double rho_line;
        for (int ai = 0; ai < size[0]; ai++) {
            tmp = 0;
            for (int ci = 0; ci < size[2]; ci++) {
                for (int bi = 0; bi < size[1]; bi++) {
                    tmp += diff_data.at(ai, bi, ci);
                }
            }
            rho_line = tmp * cell_volume_per_unit / pow(bohr_to_angstrom, 3) / len_ai;
            diff_data_red_a.push_back(rho_line);
        }
    }

    std::cout << "Preparing to output the data\n";
    // output the data and make the plot
    if (std::find(abscissa.begin(), abscissa.end(), "c") != abscissa.end()) {
        std::ofstream fout;
        fout.open(output_file+".diff.1d.c.data");
        fout << "#c(angstrom) diff_rho(e) (number of electron per Angstrom)\n";
        arma::vec c_coord = arma::linspace(0, c, diff_data_red_c.size());
        for (int i = 0; i < diff_data_red_c.size(); i++) {
            fout << c_coord.at(i) << " " << diff_data_red_c.at(i) << "\n";
        }
        fout.close();
    }
    
    if (std::find(abscissa.begin(), abscissa.end(), "b") != abscissa.end()) {
        std::ofstream fout;
        fout.open(output_file+".diff.1d.b.data");
        fout << "#b(angstrom) diff_rho(e) (number of electron per Angstrom)\n";
        arma::vec b_coord = arma::linspace(0, b, diff_data_red_b.size());
        for (int i = 0; i < diff_data_red_b.size(); i++) {
            fout << b_coord.at(i) << " " << diff_data_red_b.at(i) << "\n";
        }
        fout.close();
    }

    if (std::find(abscissa.begin(), abscissa.end(), "a") != abscissa.end())  {
        std::ofstream fout;
        fout.open(output_file+".diff.1d.a.data");
        fout << "#a(angstrom) diff_rho(e) (number of electron per Angstrom)\n";
        arma::vec a_coord = arma::linspace(0, a, diff_data_red_a.size());
        for (int i = 0; i < diff_data_red_a.size(); i++) {
            fout << a_coord.at(i) << " " << diff_data_red_a.at(i) << "\n";
        }
        fout.close();
    }

    return 1;
}


} // namespace atomsciflow
