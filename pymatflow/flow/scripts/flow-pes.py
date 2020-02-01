#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import argparse
import pymatgen as mg

from pymatflow.qe.static import static_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync


import matplotlib.pyplot as plt

"""
usage:
"""


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

    

    # ---------------------------------------------------------------
    # for PES 
    # ---------------------------------------------------------------
    parser.add_argument("--last-n-move", type=int,
            default=1,
            help="the last n atom in the structure file will move.")
    parser.add_argument("--xrange", type=float, nargs="+",
            default=[1, 3, 0.5],
            help="x range for moving the specified moving atoms.")
    parser.add_argument("--yrange", type=float, nargs="+",
            default=[3, 5, 0.5],
            help="y range for moving the specified moving atoms.")
    parser.add_argument("--zshift", type=float,
            default=0.0,
            help="z shift for the moving atoms, will shift the z of specified moving atoms by value of zshift")

    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=0,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="pes-search",
            help="jobname on the pbs server")
    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")
    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


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


    if os.path.exists(args.directory):
        shutil.rmtree(args.directory)
    os.mkdir(args.directory)

    shutil.copyfile(task.arts.xyz.file, os.path.join(args.directory, task.arts.xyz.file))
    all_upfs = [s for s in os.listdir() if s.split(".")[-1] == "UPF"]
    for element in task.arts.xyz.specie_labels:
        for upf in all_upfs:
            if upf.split(".")[0] == element:
                shutil.copyfile(upf, os.path.join(args.directory, upf))
                break
    # 

    os.chdir(args.directory)

    # shift z of the specified atoms by args.shfitz
    #for atom in task.arts.xyz.atoms[-args.last_n_move:]:
    #    atoms.z = atoms.z + args.zshift
    #task.arts.xyz.to_xyz(args.file.split(".xyz")[0]+".zshifted.xyz")
    # end shift z of the atoms
    #----------------------------------------------
    #@@@@ now we don't shift z here in python 
    #@@@@ we shift z in the bash script
    #----------------------------------------------

    with open("static-scf.in.template", 'w') as fout:
        task.control.to_in(fout)
        task.system.to_in(fout)
        task.electrons.to_in(fout)
        #task.ions.to_int(fout)

        fout.write("ATOMIC_SPECIES\n")
        upf_all = [s for s in os.listdir("./") if s.split(".")[-1] == "UPF"]
        for element in task.arts.xyz.specie_labels:
            for upf in upf_all:
                if upf.split(".")[0] == element:
                    pseudo_file =upf
                    break
            fout.write("%s %f %s\n" % (element, mg.Element(element).atomic_mass, pseudo_file))
            pseudo_file = None
            # after pseudo_file used, set it to None to avoid it will be used in the next element
        fout.write("\n")
        
        # write cell
        cell = task.arts.xyz.cell
        fout.write("CELL_PARAMETERS angstrom\n")
        for i in range(3):
            fout.write("%.9f %.9f %.9f\n" % (cell[i][0], cell[i][1], cell[i][2]))
        fout.write("\n")

        # writing KPOINTS to the fout
        task.arts.write_kpoints(fout)
        # =========================
        # writing forces act on atoms
        if task.arts.atomic_forces_status == True:
            task.arts.write_atomic_forces(fout)
        # =========================


    # write job control script
    with open("pes.pbs", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("#PBS -N %s\n" % "pes")
        fout.write("#PBS -l nodes=%d:ppn=%d\n" % (1, 32))
        fout.write("\n")
        fout.write("cd $PBS_O_WORKDIR\n")
        fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
      
        fout.write("# z of atoms to move will be shifted by zshift\n")
        fout.write("zshift=%f\n" % args.zshift)
        fout.write("\n")
        fout.write("# last_n_move atoms to move\n")
        fout.write("last_n_move=%d\n" % args.last_n_move)
        fout.write("for deltax in `seq -w %.3f %.3f %.3f`\n" % (args.xrange[0], args.xrange[2], args.xrange[1]))
        fout.write("do\n")
        fout.write("for deltay in `seq -w %.3f %.3f %.3f`\n" % (args.yrange[0], args.yrange[2], args.yrange[1]))
        fout.write("do\n")
        fout.write("  mkdir -p _${deltax}_${deltay}_\n")
        fout.write("  cp *.UPF _${deltax}_${deltay}_\n")
        fout.write("  cp static-scf.in.template _${deltax}_${deltay}_/static-scf.in\n")
        fout.write("  cat >> _${deltax}_${deltay}_/static-scf.in<<EOF\n")
        fout.write("\n")
        fout.write("ATOMIC_POSITIONS angstrom\n")
        fout.write("EOF\n")
        fout.write("  cat %s | tail -n +3 | head -n -${last_n_move} >> _${deltax}_${deltay}_/static-scf.in\n" % (args.file))
        fout.write("  # read the last n atoms and modify the x y z of them and put them into the input file\n")
        fout.write("  for i in `seq -w 1 1 ${last_n_move}`\n")
        fout.write("  do\n")
        fout.write("    lastn=`echo \"${last_n_move} + 1 - ${i}\" | bc`\n")
        fout.write("    atom=`cat %s | tail -n -${last_n_move} | tail -n -${lastn} | head -n 1`\n" % (args.file))
        fout.write("    label=`echo $atom | cut -d ' ' -f1`\n")
        fout.write("    x=`echo $atom | cut -d ' ' -f2`\n")
        fout.write("    y=`echo $atom | cut -d ' ' -f3`\n")
        fout.write("    z=`echo $atom | cut -d ' ' -f4`\n")
        fout.write("    x=`echo \"${x} + ${deltax}\" | bc`\n")
        fout.write("    y=`echo \"${y} + ${deltay}\" | bc`\n")
        fout.write("    z=`echo \"${z} + ${zshift}\" | bc`\n")
        fout.write("    cat >> _${deltax}_${deltay}_/static-scf.in<<EOF\n")
        fout.write("${label} ${x} ${y} ${z}\n")
        fout.write("EOF\n")
        fout.write("  done\n")
        fout.write("  # run the calculation\n")
        fout.write("  cd _${deltax}_${deltay}_\n")
        fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % ("pw.x", "static-scf.in", "static-scf.out"))
        fout.write("  cd ../\n")
        fout.write("done\n")
        fout.write("done\n")
    
    # write result analysis file
    with open("get_pes.sh", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("cat > pes.data<<EOF\n")
        fout.write("# format: x y energy(Ry)\n")
        fout.write("EOF\n")
        fout.write("\n")
        fout.write("for deltax in `seq -w %.3f %.3f %.3f`\n" % (args.xrange[0], args.xrange[2], args.xrange[1]))
        fout.write("do\n")
        fout.write("for deltay in `seq -w %.3f %.3f %.3f`\n" % (args.yrange[0], args.yrange[2], args.yrange[1]))
        fout.write("do\n")
        fout.write("  energy=`cat _${deltax}_${deltay}_/static-scf.out | grep '!    total energy' | tail -1`\n")
        fout.write("  cat >> pes.data<<EOF\n")
        fout.write("${deltax} ${deltay} ${energy:32:-2}\n")
        fout.write("EOF\n")
        fout.write("done\n")
        fout.write("done\n")
        fout.write("\n")
        fout.write("cat > plot.gnuplot<<EOF\n")
        fout.write("set term png\n")
        fout.write("set output 'pes.png'\n")
        fout.write("splot 'pes.data'\n")
        fout.write("EOF\n")
        fout.write("gnuplot plot.gnuplot\n")
    
    os.chdir("../")



    # server handle
    if args.auto == 0:
        pass
    elif args.auto == 1:
        mover = rsync()
        if args.server == "pbs":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_pbs.conf"))
        elif args.server == "yh":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_yh.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
    elif args.auto == 2:
        mover = rsync()
        if args.server == "pbs":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_pbs.conf"))
        elif args.server == "yh":
            mover.get_info(os.path.join(os.path.expanduser("~"), ".pymatflow/server_yh.conf"))
        mover.copy_default(source=os.path.abspath(args.directory))
        ctl = ssh()
        if args.server == "pbs":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_pbs.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="pes.pbs", server="pbs")
        elif args.server == "yh":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="pes.sub", server="yh")
