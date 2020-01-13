#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import argparse

from pymatflow.qe.static import static_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync


import matplotlib.pyplot as plt

"""
usage:
    qe-scf.py -f xxx.xyz -k '2 2 2 0 0 0' --ecutwfc 100
"""


def pes(task, directory, runopt, mpi):
    stridex = [1.0, 5.0]
    stridey = [3.0, 5.0]
    stepx = 0.5
    stepy = 0.5
    sites = []
    for i in range(int((stridex[1]-stridex[0])/stepx)):
        for j in  range(int((stridey[1]-stridey[0])/stepy)):
            #sites.append([stridex[0]+i*stepx, stridey[0]+j*stepy])
            sites.append({"xy": [stridex[0]+i*stepx, stridey[0]+j*stepy], "energy": None})

#    if os.path.exists(directory):
#        shutil.rmtree(directory)
#    os.mkdir(directory)
#    os.system("cp *.UPF ./%s/" % directory) 
    os.chdir(directory)

    for site in sites:
#        task.arts.xyz.atoms[-1].x = site["xy"][0]
#        task.arts.xyz.atoms[-1].y = site["xy"][1]
#        task.arts.xyz.to_xyz("lih-slab-h.xyz")
        
#        task.scf(directory="%.3f-%.3f" % (site["xy"][0], site["xy"][1]), runopt=runopt, mpi=mpi)

        with open("%.3f-%.3f/static-scf.out" % (site["xy"][0], site["xy"][1]), 'r') as fin:
            lines = fin.readlines()

        for line in lines:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "!" and line.split()[5] == "Ry":
                site["energy"] = float(line.split()[4])
    ax = plt.axes(projection='3d')
    x = [sites[i]["xy"][0] for i in range(len(sites))]
    y = [sites[i]["xy"][1] for i in range(len(sites))]
    energy = [sites[i]["energy"] for i in range(len(sites))]
    with open("pes.data", 'w') as fout:
        for i in range(len(sites)):
            fout.write("%.9f %.9f %.9f\n" % (sites[i]["xy"][0], sites[i]["xy"][1], sites[i]["energy"]))
    ax.scatter(x, y, energy)
    plt.savefig("pes.matplotlib.png")
    plt.show()

control = {}
system = {}
electrons = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-pes-static",
            help="Directory for the static running.")
    parser.add_argument("-f", "--file", type=str,
            help="The xyz file name.")
    parser.add_argument("--runopt", type=str, default="genrun", 
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")
    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    # -------------------------------------------------------------------
    #                       scf related parameters
    # -------------------------------------------------------------------
    parser.add_argument("--ecutwfc", type=int, default=100,
            help="Kinetic energy cutoff for wave functions in unit of Rydberg, default value: 100 Ry")

    parser.add_argument("--ecutrho", type=int, default=400,
            help="Kinetic energy cutoff for charge density and potential in unit of Rydberg, default value: 400 Ry")

    parser.add_argument("--kpoints-option", type=str, default="automatic", 
            choices=["automatic", "gamma", "tpiba_b"],
            help="Kpoints generation scheme option for the SCF or non-SCF calculation")

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 1 1 1 0 0 0")

    parser.add_argument("--conv-thr", type=float, default=1.0e-6,
            help="Convergence threshold for SCF calculation.")

    parser.add_argument("--occupations", type=str, default="smearing",
            choices=["smearing", "tetrahedra", "tetrahedra_lin", "tetrahedra_opt", "fixed", "from_input"],
            help="Occupation method for the calculation.")
    
    parser.add_argument("--smearing", type=str, default="gaussian",
            choices=["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"],
            help="Smearing type for occupations by smearing, default is gaussian in this script")

    parser.add_argument("--degauss", type=float, default=0.001,
            help="Value of the gaussian spreading (Ry) for brillouin-zone integration in metals.(defualt: 0.001 Ry)")

    parser.add_argument("--vdw-corr", type=str, default="none",
            choices=["dft-d", "dft-d3", "ts", "xdm"],
            help="Type of Van der Waals correction in the calculation")

    parser.add_argument("--nbnd", type=int, default=None,
            help="Number of electronic states (bands) to be calculated")

    parser.add_argument("--tstress", type=str, default=".false.",
            choices=[".true.", ".false."],
            help="calculate stress. default=.false.")

    # -----------------------------------------------------------
    #           ATOMIC_FORCES
    # -----------------------------------------------------------
    parser.add_argument("--pressure", type=float, default=None,
            help="specify pressure acting on system in unit of Pa")
    parser.add_argument("--pressuredir", type=str, default=None,
            choices=["x", "y", "z"],
            help="specify direction of pressure acting on system.")
    
    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.emuhelper/server.conf")


    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================   
    args = parser.parse_args()
    xyzfile = args.file
    control["tstress"] = args.tstress
    system["ecutwfc"] = args.ecutwfc
    system["ecutrho"] = args.ecutrho
    system["occupations"] = args.occupations
    system["smearing"] = args.smearing
    system["degauss"] = args.degauss
    system["vdw_corr"] = args.vdw_corr
    system["nbnd"] = args.nbnd
    electrons["conv_thr"] = args.conv_thr
    

    task = static_run()
    task.get_xyz(xyzfile)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
    task.set_params(control=control, system=system, electrons=electrons)
    #task.set_atomic_forces(pressure=args.pressure, pressuredir=args.pressuredir)
    pes(task, directory=args.directory, runopt=args.runopt, mpi=args.mpi)
    
    #task.scf(directory=args.directory, runopt=args.runopt, mpi=args.mpi)

