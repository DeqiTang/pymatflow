#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import argparse

import pymatflow.base as base

from pymatflow.flow.surface_pes import qe_run
#from pymatflow.remote.server import server_handle


"""
usage:
"""


control = {}
system = {}
electrons = {}
ions = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-pes-relax",
            help="Directory for the pes relax running.")

    parser.add_argument("--xyz", type=str,
            help="The xyz file name.")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")


    # -------------------------------------------------------------------
    #                       scf related parameters
    # -------------------------------------------------------------------
    parser.add_argument("--ecutwfc",
            type=int, default=100)

    parser.add_argument("--ecutrho",
            type=int, default=400)

    parser.add_argument("--kpoints-option", type=str, default="automatic",
            choices=["automatic", "gamma", "crystal_b"],
            help="Kpoints generation scheme option for the SCF or non-SCF calculation")

    parser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 1 1 1 0 0 0")

    parser.add_argument("--conv-thr", type=float, default=1.0e-6,
            help="the conv_thr for scf, when doing geometric optimization better use a strict covnergec for scf")

    parser.add_argument("--occupations", type=str, default="smearing",
            choices=["smearing", "tetrahedra", "tetrahedra_lin", "tetrahedra_opt", "fixed", "from_input"],
            help="Occupation method for the calculation.")

    parser.add_argument("--smearing", type=str, default="gaussian",
            choices=["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"],
            help="Smearing type for occupations by smearing, default is gaussian in this script")

    parser.add_argument("--degauss", type=float, default=0.001,
            help="Value of the gaussian spreading (Ry) for brillouin-zone integration in metals.(defualt: 0.001 Ry)")

    parser.add_argument("--vdw-corr", help="vdw_corr = dft-d, dft-d3, ts-vdw, xdm", type=str, default="none")


    # -------------------------------------------------------------------
    #               geometric optimization related parameters
    # -------------------------------------------------------------------
    parser.add_argument("--etot-conv-thr",
            type=float, default=1.0e-4,
            help="convergence threshold of energy for geometric optimization")

    parser.add_argument("--forc-conv-thr",
            type=float, default=1.0e-3,
            help="convergence threshold for force in optimization,(usually it is more important than energy)")

    parser.add_argument("--nstep",
            type=int, default=50,
            help="maximum ion steps for geometric optimization")


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
    parser.add_argument("--move-atom", type=int, nargs="+",
            default=[-1],
            help="specify the atoms to move.")

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
    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="pes-relax",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")

    parser.add_argument("--queue", type=str, default=None,
            help="queue to submit the job")


    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    control["etot_conv_thr"] = args.etot_conv_thr
    control["forc_conv_thr"] = args.forc_conv_thr
    control["nstep"] = args.nstep
    system["ecutwfc"] = args.ecutwfc
    system["ecutrho"] = args.ecutrho
    system["occupations"] = args.occupations
    system["smearing"] = args.smearing
    system["degauss"] = args.degauss
    system["vdw_corr"] = args.vdw_corr
    electrons["conv_thr"] = args.conv_thr

    task = qe_run()
    task.get_xyz(args.xyz)
    task.set_relax()
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
    task.set_params(control=control, system=system, electrons=electrons, ions=ions)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn, queue=args.queue)
    task.set_pes(move_atom=args.move_atom, xrange=args.xrange, yrange=args.yrange, zshift=args.zshift)
    task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)


"""
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

    self.arts.pseudo.dir = os.path.abspath(directory)
    self.control.pseudo_dir = os.path.abspath(directory)
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
            fout.write("%s %f %s\n" % (element, base.element[element].mass, pseudo_file))
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
    server_handle(auto=args.auto, directory=args.directory, jobfilebase="pes-relax", server=args.server)
"""
