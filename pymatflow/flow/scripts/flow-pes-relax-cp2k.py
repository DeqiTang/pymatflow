#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import argparse
import pymatgen as mg

from pymatflow.cp2k.opt import opt_run
from pymatflow.remote.ssh import ssh
from pymatflow.remote.rsync import rsync


import matplotlib.pyplot as plt

"""
usage:
"""



params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="directory of the calculation", type=str, default="tmp-cp2k-pes-opt")
    parser.add_argument("-f", "--file", help="the xyz file name", type=str)

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--runopt", type=str, default="genrun", 
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("-p", "--printout-option", nargs="+", type=int,
            default=[],
            help="Properties printout option(0, 1, 2 implemented now), you can also activate multiple prinout-option at the same time")

    parser.add_argument("--ls-scf", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="use linear scaling scf method")
    
    parser.add_argument("--qs-method", type=str, default="GPW",
            choices=["AM1", "DFTB", "GAPW", "GAPW_XC", "GPW", "LRIGPW", "MNDO", "MNDOD", 
                "OFGPW", "PDG", "PM3", "PM6", "PM6-FM", "PNNL", "RIGPW", "RM1"],
            help="specify the electronic structure method")

    parser.add_argument("--eps-scf", type=float, default=1.0e-6,
            help="DFT-SCF-EPS_SCF")

    parser.add_argument("--xc-functional", type=str, default="PBE",
            help="XC_FUNCTIONAL type")

    parser.add_argument("--cutoff", type=int, default=100,
            help="CUTOFF, default value: 100 Ry, if you find your SCF hard to converge, you can try increasing CUTOFF")

    parser.add_argument("--rel-cutoff", type=int, default=60,
            help="REL_CUTOFF, default value: 60 Ry")

    parser.add_argument("-k", "--kpoints-scheme", type=str,
            default="GAMMA",
            help="DFT-KPOINTS-SCHEME(str): can be NONE, GAMMA, MONKHORST-PACK, MACDONALD, GENERAL. when you set MONKHORST-PACK, you should also add the three integers like 'monkhorst-pack 3 3 3'")
    
    parser.add_argument("--diag", type=str, default="TRUE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing tranditional diagonalization for SCF")

    parser.add_argument("--ot", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing orbital transformation for SCF")

    parser.add_argument("--alpha", type=float, default=0.4,
            help="DFT-SCF-MIXING-ALPHA")

    parser.add_argument("--smear", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="switch on or off smearing for occupation")

    parser.add_argument("--added-mos", help="ADDED_MOS for SCF", type=int, default=0)

    parser.add_argument("--smear-method", help="smearing type: FERMI_DIRAC, ENERGY_WINDOW", type=str, default="FERMI_DIRAC")

    parser.add_argument("--electronic-temp", help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR", type=float, default=300)

    parser.add_argument("--window-size", help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing", type=float, default=0)

    # motion/geo_opt related
    parser.add_argument("--optimizer", type=str, default="BFGS",
            help="optimization algorithm for geometry optimization: BFGS, CG, LBFGS")
    parser.add_argument("--max-iter", type=int, default=200,
            help="maximum number of geometry optimization steps.")
    parser.add_argument("--type", type=str, default="MINIMIZATION",
            help="specify which kind of geometry optimization to perform: MINIMIZATION(default), TRANSITION_STATE")
    parser.add_argument("--max-dr", type=float, default=3e-3,
            help="Convergence criterion for the maximum geometry change between the current and the last optimizer iteration.")
    parser.add_argument("--max-force", type=float, default=4.50000000E-004,
            help="Convergence criterion for the maximum force component of the current configuration.")
    parser.add_argument("--rms-dr", type=float, default=1.50000000E-003,
            help="Convergence criterion for the root mean square (RMS) geometry change between the current and the last optimizer iteration.")
    parser.add_argument("--rms-force", type=float, default=3.00000000E-004,
            help="Convergence criterion for the root mean square (RMS) force of the current configuration.")


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
    parser.add_argument("--jobname", type=str, default="pes-opt",
            help="jobname on the pbs server")
    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")
    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")




    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================   
    args = parser.parse_args()
    
    params["FORCE_EVAL-DFT-LS_SCF"] = args.ls_scf
    params["FORCE_EVAL-DFT-QS-METHOD"] = args.qs_method
    params["FORCE_EVAL-DFT-MGRID-CUTOFF"] = args.cutoff
    params["FORCE_EVAL-DFT-MGRID-REL_CUTOFF"] = args.rel_cutoff
    params["FORCE_EVAL-DFT-XC-XC_FUNCTIONAL"] = args.xc_functional
    params["FORCE_EVAL-DFT-SCF-EPS_SCF"] = args.eps_scf
    params["FORCE_EVAL-DFT-SCF-ADDED_MOS"] = args.added_mos
    params["FORCE_EVAL-DFT-SCF-SMEAR"] = args.smear
    params["FORCE_EVAL-DFT-SCF-SMEAR-METHOD"] = args.smear_method
    params["FORCE_EVAL-DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"] = args.electronic_temp
    params["FORCE_EVAL-DFT-SCF-SMEAR-WINDOW_SIZE"] = args.window_size
    params["FORCE_EVAL-DFT-SCF-DIAGONALIZATION"] = args.diag
    params["FORCE_EVAL-DFT-SCF-OT"] = args.ot
    params["FORCE_EVAL-DFT-SCF-MIXING-ALPHA"] = args.alpha
    params["FORCE_EVAL-DFT-KPOINTS-SCHEME"] = args.kpoints_scheme

    params["MOTION-GEO_OPT-MAX_ITER"] = args.max_iter
    params["MOTION-GEO_OPT-OPTIMIZER"] = args.optimizer
    params["MOTION-GEO_OPT-TYPE"] = args.type
    params["MOTION-GEO_OPT-MAX_DR"] = args.max_dr
    params["MOTION-GEO_OPT-MAX_FORCE"] = args.max_force
    params["MOTION-GEO_OPT-RMS_DR"] = args.rms_dr
    params["MOTION-GEO_OPT-RMS_FORCE"] = args.rms_force

    task = opt_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    #task.geo_opt(directory=args.directory, mpi=args.mpi, runopt=args.runopt)
    

    if os.path.exists(args.directory):
        shutil.rmtree(args.directory)
    os.mkdir(args.directory)

    shutil.copyfile(task.force_eval.subsys.xyz.file, os.path.join(args.directory, task.force_eval.subsys.xyz.file))
    # 

    os.chdir(args.directory)

    # shift z of the specified atoms by args.shfitz
    #----------------------------------------------
    #@@@@ now we don't shift z here in python 
    #@@@@ we shift z in the bash script
    #----------------------------------------------

    with open("relax.in.template", 'w') as fout:
        task.control.to_in(fout)
        task.system.to_in(fout)
        task.electrons.to_in(fout)
        task.ions.to_in(fout)

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
    with open("pes-relax.pbs", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("#PBS -N %s\n" % "pes")
        fout.write("#PBS -l nodes=%d:ppn=%d\n" % (args.nodes, args.ppn))
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
        fout.write("  cp relax.in.template _${deltax}_${deltay}_/relax.in\n")
        fout.write("  cat >> _${deltax}_${deltay}_/relax.in<<EOF\n")
        fout.write("\n")
        fout.write("ATOMIC_POSITIONS angstrom\n")
        fout.write("EOF\n")
        #fout.write("  cat %s | tail -n +3 | head -n -${last_n_move} >> _${deltax}_${deltay}_/relax.in\n" % (args.file))
        fout.write("  # read the main structure atoms(unmoved) and put them into the input file\n")
        fout.write("  natom=`cat %s | head -n 1`\n" % args.file)
        fout.write("  end_of_fix=`echo \"${natom} - ${last_n_move} + 2\" | bc`\n")
        fout.write("  for i in `seq -w 3 1 ${end_of_fix}`\n")
        fout.write("  do\n")
        fout.write("    atom=`cat %s | head -n ${i} | tail -n 1`\n" % (args.file))
        fout.write("    label=`echo $atom | cut -d ' ' -f1`\n")
        fout.write("    x=`echo $atom | cut -d ' ' -f2`\n")
        fout.write("    y=`echo $atom | cut -d ' ' -f3`\n")
        fout.write("    z=`echo $atom | cut -d ' ' -f4`\n")
        fout.write("    cat >> _${deltax}_${deltay}_/relax.in<<EOF\n")
        fout.write("${label} ${x} ${y} ${z} 0 0 0\n")
        fout.write("EOF\n")
        fout.write("  done\n")
        fout.write("\n")
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
        fout.write("    cat >> _${deltax}_${deltay}_/relax.in<<EOF\n")
        fout.write("${label} ${x} ${y} ${z} 0 0 1\n")
        fout.write("EOF\n")
        fout.write("  done\n")
        fout.write("  # run the calculation\n")
        fout.write("  cd _${deltax}_${deltay}_\n")
        fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % ("pw.x", "relax.in", "relax.out"))
        fout.write("  cd ../\n")
        fout.write("done\n")
        fout.write("done\n")


    # write bash script to generate the xyz trajectory file -> (unrelaxed original traj)
    # note: bash script runs slowly to generate the trajectory file
    #       we might better use python to directly genrate the trajectory file
    with open("get_traj_unrelaxed.sh", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("\n")
        fout.write("output_trajfile=trajectory-unrelaxed.xyz\n")
        fout.write("# z of atoms to move will be shifted by zshift\n")
        fout.write("zshift=%f\n" % args.zshift)
        fout.write("\n")
        fout.write("# last_n_move atoms to move\n")
        fout.write("last_n_move=%d\n" % args.last_n_move)
        fout.write("for deltax in `seq -w %.3f %.3f %.3f`\n" % (args.xrange[0], args.xrange[2], args.xrange[1]))
        fout.write("do\n")
        fout.write("for deltay in `seq -w %.3f %.3f %.3f`\n" % (args.yrange[0], args.yrange[2], args.yrange[1]))
        fout.write("do\n")
        fout.write("  cat %s | head -n 1 >> ${output_trajfile}\n" % args.file)
        fout.write("  cat >> ${output_trajfile}<<EOF\n")
        fout.write("${deltax} ${deltay}\n")
        fout.write("EOF\n")
        fout.write("  cat %s | tail -n +3 | head -n -${last_n_move} >> ${output_trajfile}\n" % (args.file))
        fout.write("  # read the last n atoms and modify the x y z of them and put them into the trajectory file\n")
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
        fout.write("    cat >> ${output_trajfile}<<EOF\n")
        fout.write("${label} ${x} ${y} ${z}\n")
        fout.write("EOF\n")
        fout.write("  done\n")
        fout.write("done\n")
        fout.write("done\n")

    # write bash script to generate the xyz trajectory file -> (relaxed)
    with open("get_traj_relaxed.sh", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("\n") 
        fout.write("\n")
        fout.write("output_trajfile=trajectory-relaxed.xyz\n")
        fout.write("natom=`cat %s | head -n 1`\n" % args.file)
        fout.write("for deltax in `seq -w %.3f %.3f %.3f`\n" % (args.xrange[0], args.xrange[2], args.xrange[1]))
        fout.write("do\n")
        fout.write("for deltay in `seq -w %.3f %.3f %.3f`\n" % (args.yrange[0], args.yrange[2], args.yrange[1]))
        fout.write("do\n")
        fout.write("  echo ${natom} >> ${output_trajfile}\n")
        fout.write("  cat >> ${output_trajfile}<<EOF\n")
        fout.write("${deltax} ${deltay}\n")
        fout.write("EOF\n")
        fout.write("  begin_final=`cat _${deltax}_${deltay}_/relax.out | grep -n \"Begin final\" | cut -d ':' -f1`\n")
        fout.write("  end_final=`cat _${deltax}_${deltay}_/relax.out | grep -n \"End final\" | cut -d ':' -f1`\n")
        fout.write("  cat _${deltax}_${deltay}_/relax.out | head -n +${end_final} | tail -n +${begin_final} | head -n -1 | tail -n +4 >> ${output_trajfile}\n")
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
        fout.write("  energy=`cat _${deltax}_${deltay}_/relax.out | grep '!    total energy' | tail -1`\n")
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
            ctl.submit(workdir=args.directory, jobfile="pes-relax.pbs", server="pbs")
        elif args.server == "yh":
            ctl.get_info(os.path.join(os.path.expanduser('~'), ".pymatflow/server_yh.conf"))
            ctl.login()
            ctl.submit(workdir=args.directory, jobfile="pes-relax.sub", server="yh")
