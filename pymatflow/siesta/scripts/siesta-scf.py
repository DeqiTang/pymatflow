#!/usr/bin/evn python
# _*_ coding:utf-8 _*_

import argparse

from pymatflow.siesta.static import static_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync

"""
usage:
    qe-scf.py -f xxx.xyz -k '2 2 2 0 0 0' -- 100
"""


electrons = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-siesta-static")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)

    parser.add_argument("--runopt", type=str, default="genrun", 
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--mpi", help="MPI command", type=str, default="")
    parser.add_argument("--meshcutoff", help="MeshCutoff (Ry)", type=int, default=200)
    parser.add_argument("--solution-method", help="SolutionMethod(diagon, OMM, OrderN, PEXSI)", type=str, default="diagon")
    parser.add_argument("--functional", help="XC.functional", type=str, default="GGA")
    parser.add_argument("--authors", help="XC.authors", type=str, default="PBE")
    parser.add_argument("--tolerance", help="DM.Tolerance", type=float, default=1.0e-6)
    parser.add_argument("--numberpulay", help="DM.NumberPulay", type=int ,default=8)
    parser.add_argument("--mixing", help="DM.MixingWeight", type=float, default=0.1)

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[3, 3, 3],
            help="set kpoints like '3 3 3'")

    parser.add_argument("--occupation", help="OccupationFunction(FD or MP)", type=str, default="FD")
    parser.add_argument("--electronic-temperature", help="Electronic Temperature", type=int, default=300)
    
    
    # ------------------------------
    # properties related parameter
    # ------------------------------
    parser.add_argument("-p", "--properties" ,nargs="+", type=int, default=[],
            help="Option for properties calculation")

    parser.add_argument("--pdos-block", type=float, nargs="+",
            default=[-20, 10, 0.2, 500])
    #------------------------------------------------------------------------------------------------
    #parser.add_argument("--bandlines", nargs="+", type=str,
    #        default=["1 0.0 0.0 0.0 \Gamma", "20 1.0 1.0 1.0 L", "20 2.0 0.0 0.0 X"],
    #        help="BandLines for band structre calculation(either choose BandLines or BandPoints)")
    #parser.add_argument("--bandpoints", nargs="+", type=str,
    #        default=["0.0 0.0 0.0", "1.0 0.0 0.0", "0.5 0.5 0.5"],
    #        help="BandPoints for band structure calculation(either choose BandPoints or BandLines)")
    # we now use seekapth to calculate the BandLines automatically so it is not needed to be set by 
    # user manually.
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

    # for server
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")

    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    xyzfile = args.file
    directory = args.directory
    
    electrons["MeshCutoff"] = args.meshcutoff
    electrons["SolutionMethod"] = args.solution_method
    electrons["XC.funtional"] = args.functional
    electrons["XC.authors"] = args.authors
    electrons["DM.Tolerance"] = args.tolerance
    electrons["DM.NumberPulay"] = args.numberpulay
    electrons["DM.MixingWeight"] = args.mixing
    electrons["OccupationFunction"] = args.occupation
    electrons["ElectronicTemperature"] = args.electronic_temperature


    task = static_run()
    task.get_xyz(xyzfile)

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

    task.scf(directory=directory, runopt=args.runopt, mpi=args.mpi, electrons=electrons, properties=args.properties, kpoints_mp=args.kpoints_mp)

    # server handle
    if args.auto == 0:
        pass
    elif args.auto == 1:
        mover = rsync()
        mover.get_info(os.path.join(os.path.expanduser("~"), ".emuhelper/server.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
    elif args.auto == 2:
        mover = rsync()
        mover.get_info(os.path.join(os.path.expanduser("~"), ".emuhelper/server.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
        ctl = ssh()
        ctl.get_info(os.path.join(os.path.expanduser('~'), ".emuhelper/server.conf"))
        ctl.login()
        ctl.submit(workdir=args.directory, jobfile="static-scf.fdf.sub")
