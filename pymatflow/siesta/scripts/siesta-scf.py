#!/usr/bin/evn python
# _*_ coding:utf-8 _*_

import os
import argparse

from pymatflow.siesta.static import static_run

"""
usage:
    qe-scf.py -f xxx.xyz -k '2 2 2 0 0 0' -- 100
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-siesta-static",
            help="directory to do the static scf calculation")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz structre file with second line specifying cell parameters")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    # --------------------------------------------------------------------------
    parser.add_argument("--meshcutoff", type=int, default=200,
            help="MeshCutoff (Ry)")

    parser.add_argument("--solution-method", type=str, default="diagon",
            choices=["diagon", "OMM", "OrderN", "PEXSI"],
            help="SolutionMethod(diagon, OMM, OrderN, PEXSI)")

    parser.add_argument("--functional", type=str, default="GGA",
            help="XC.functional")

    parser.add_argument("--authors", type=str, default="PBE",
            help="XC.authors")

    parser.add_argument("--tolerance", type=float, default=1.0e-6,
            help="DM.Tolerance")

    parser.add_argument("--numberpulay", type=int, default=8,
            help="DM.NumberPulay")

    parser.add_argument("--mixing", type=float, default=0.1,
            help="DM.MixingWeight")

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[3, 3, 3],
            help="set kpoints like '3 3 3'")

    parser.add_argument("--occupation", type=str, default="FD",
            choices=["FD", "MP"],
            help="OccupationFunction(FD or MP)")

    parser.add_argument("--electronic-temperature", type=int, default=300,
            help="Electronic Temperature")


    # ------------------------------
    # properties related parameter
    # ------------------------------
    parser.add_argument("-p", "--properties" ,nargs="+", type=int, default=[],
            help="Option for properties calculation")

    parser.add_argument("--pdos-block", type=float, nargs="+",
            default=[-20, 10, 0.2, 500])
    #------------------------------------------------------------------------------------------------
    parser.add_argument("--bandlines", nargs="+", type=str,
            help="BandLines for band structre calculation")

    parser.add_argument("--bandlines-file", type=str, default="kpath-from-seekpath.txt",
            help="BandLines for band structure calculation from file")

    #------------------------------------------------------------------------------------------------
    parser.add_argument("--polarization-grids", nargs="+", type=str,
            default=["10 3 3 no", "2 20 2 no", "4 4 15 no"],
            help="PolarizationGrids")

    parser.add_argument("--external-electric-field", nargs="+", type=float,
            default=[0.0, 0.0, 0.5],
            help="External Electric field")

    parser.add_argument("--optical-energy-minimum", type=float,
            default=0.0,
            help="Optical.Energy.Minimum")

    parser.add_argument("--optical-energy-maximum", type=float,
            default=10.0,
            help="Optical.Energy.Maximum")

    parser.add_argument("--optical-broaden", type=float,
            default=0.0,
            help="Optical.Broaden")

    parser.add_argument("--optical-scissor", type=float,
            default=0.0,
            help="Optical.Scissor")

    parser.add_argument("--optical-mesh", nargs="+", type=int,
            default=[5, 5, 5],
            help="Optical.Mesh")

    parser.add_argument("--optical-polarization-type", type=str,
            default="unpolarized",
            help="Optical.PolarizationType")

    parser.add_argument("--optical-vector", nargs="+", type=float,
            default=[1.0, 0.0, 0.5],
            help="Optical.Vector")

    parser.add_argument("--wannier90-unkgrid", nargs="+", type=int,
            default=[10, 10, 10],
            help="Siesta2Wannier90.UnkGrid[1-3]")


    # -----------------------------------------------------------------
    #                      run param
    # -----------------------------------------------------------------
    parser.add_argument("--mpi", type=str, default="",
            help="MPI command")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="siesta-scf",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    params = {}

    params["MeshCutoff"] = args.meshcutoff
    params["SolutionMethod"] = args.solution_method
    params["XC.funtional"] = args.functional
    params["XC.authors"] = args.authors
    params["DM.Tolerance"] = args.tolerance
    params["DM.NumberPulay"] = args.numberpulay
    params["DM.MixingWeight"] = args.mixing
    params["OccupationFunction"] = args.occupation
    params["ElectronicTemperature"] = args.electronic_temperature

    # if band structure is in the properties get the bandlines
    if 3 in args.properties and args.bandlines != None:
        # bandlines from script argument args.bandlines
        bandlines = []
        for kpoint in args.bandlines:
            if kpoint.split()[4] != "|":
                bandlines.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    int(kpoint.split()[4]),
                    ])
            elif kpoint.split()[4] == "|":
                bandlines.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    "|",
                    ])
    elif 3 in args.properties and args.bandlines == None:
        # bandlines read from file specified by args.bandlines_file
        # file is in format like this
        """
        5
        0.0 0.0 0.0 #GAMMA 15
        x.x x.x x.x #XXX |
        x.x x.x x.x #XXX 10
        x.x x.x x.x #XXX 15
        x.x x.x x.x #XXX 20
        """
        # if there is a '|' behind the label it means the path is
        # broken after that point!!!
        bandlines = []
        with open(args.bandlines_file, 'r') as fin:
            bandlines_file = fin.readlines()
        nk = int(bandlines_file[0])
        for i in range(nk):
            if bandlines_file[i+1].split("\n")[0].split()[4] != "|":
                bandlines.append([
                    float(bandlines_file[i+1].split()[0]),
                    float(bandlines_file[i+1].split()[1]),
                    float(bandlines_file[i+1].split()[2]),
                    bandlines_file[i+1].split()[3].split("#")[1].upper(),
                    int(bandlines_file[i+1].split()[4]),
                    ])
            elif bandlines_file[i+1].split("\n")[0].split()[4] == "|":
                bandlines.append([
                    float(bandlines_file[i+1].split()[0]),
                    float(bandlines_file[i+1].split()[1]),
                    float(bandlines_file[i+1].split()[2]),
                    bandlines_file[i+1].split()[3].split("#")[1].upper(),
                    '|',
                    ])
    else:
        pass
        # 3 not in args.properties
        # do not calculate the band structure
        # no need to set the bandlines
    # end the setting of bandlines

    task = static_run()
    task.get_xyz(args.file)

    task.properties.set_params(
        #bandlines = args.bandlines,
        #bandpoints = args.bandpoints,
        polarization_grids = args.polarization_grids,
        external_electric_field = args.external_electric_field,
        optical_energy_minimum = args.optical_energy_minimum,
        optical_energy_maximum = args.optical_energy_maximum,
        optical_broaden = args.optical_broaden,
        optical_scissor = args.optical_scissor,
        optical_mesh = args.optical_mesh,
        optical_polarization_type = args.optical_polarization_type,
        optical_vector = args.optical_vector,
        wannier90_unkgrid = args.wannier90_unkgrid,
        )

    if 3 in args.properties:
        task.properties.bandlines = bandlines

    task.set_params(params=params)
    task.set_kpoints(kpoints_mp=args.kpoints_mp)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto, properties=args.properties)
