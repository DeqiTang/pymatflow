#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import os
import shutil
import argparse

import pymatflow.base as base

from pymatflow.qe.opt import opt_run
from pymatflow.remote.server import server_handle

"""
"""

control = {}
system = {}
electrons = {}
ions = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory",
            type=str, default="tmp-qe-relax-cubic",
            help="directory for the relax running")

    parser.add_argument("-f", "--file",
            type=str,
            help="the xyz file containg the structure to be simulated")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--mpi",
            type=str, default="",
            help="the mpi command used")

    # -------------------------------------------------------------------
    #                       scf related parameters
    # -------------------------------------------------------------------
    parser.add_argument("--ecutwfc",
            type=int, default=100)

    parser.add_argument("--ecutrho", type=int, default=None,
            help="Kinetic energy cutoff for charge density and potential in unit of Rydberg, default value: None")

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

    # --------------------------------------------------------------------------
    # na stepa
    # --------------------------------------------------------------------------
    parser.add_argument("--na", type=int, default=10,
            help="number of a used")
    parser.add_argument("--stepa", type=float, default=0.05,
            help="a step")

    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="relax-cubic",
            help="jobname on the pbs server")
    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")
    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")



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

    task = opt_run()
    task.set_relax()
    task.get_xyz(args.file)
    task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
    task.set_params(control=control, system=system, electrons=electrons, ions=ions)
    #task.relax(directory=args.directory, runopt=args.runopt, mpi=args.mpi, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)


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

    with open("relax.in.template", 'w') as fout:
        task.control.to_in(fout)
        task.system.to_in(fout)
        task.electrons.to_in(fout)
        task.ions.to_in(fout)

        coordtype = "angstrom"
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
        if coordtype == "angstrom":
            fout.write("ATOMIC_POSITIONS angstrom\n")
            if task.arts.ifstatic == True:
                for atom in task.arts.xyz.atoms:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
            elif task.arts.ifstatic == False:
                for atom in task.arts.xyz.atoms:
                    fout.write("%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z))
                    for fix in atom.fix:
                        if fix == True:
                            fout.write("\t0")
                        elif fix == False:
                            fout.write("\t1")
                    fout.write("\n")
            else:
                print("===============================================\n")
                print("warning: qe.base.arts.to_in():\n")
                print("arts.ifstatic could only be True or False\n")
                sys.exit(1)
            fout.write("\n")
        elif coordtype == "crystal":
            # crystal namely fractional coordinate can be convert from cartesian coordinates
            # the conversion process is like transformation of presentation in quantum mechanics
            # the convmat is bulid to do the conversion
            #latcell = np.array(self.xyz.cell)
            #latcell = latcell.reshape(3, 3)
            latcell = np.array(task.arts.xyz.cell)
            convmat = np.linalg.inv(latcell.T)
            crystal_coord = np.zeros([task.arts.xyz.natom, 3])
            for i in range(task.arts.xyz.natom):
                crystal_coord[i] = convmat.dot(np.array([task.arts.xyz.atoms[i].x, task.arts.xyz.atoms[i].y, task.arts.xyz.atoms[i].z]))
            #
            fout.write("ATOMIC_POSITIONS crystal\n")
            if task.arts.ifstatic == True:
                for k in range(task.arts.xyz.natom):
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (task.arts.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
            elif task.arts.ifstatic == False:
                for k in range(task.arts.xyz.natom):
                    fout.write("%s\t%.9f\t%.9f\t%.9f" % (task.arts.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                    for fix in task.arts.xyz.atoms[k].fix:
                        if fix == True:
                            fout.write("\t0")
                        elif fix == False:
                            fout.write("\t1")
                    fout.write("\n")
            else:
                print("===============================================\n")
                print("warning: qe.base.arts.to_in():\n")
                print("arts.ifstatic could only be True or False\n")
                sys.exit(1)
            fout.write("\n")
        # end crystal type ATOMIC_POSITIONS

        # writing KPOINTS to the fout
        task.arts.write_kpoints(fout)
        # =========================
        #
        # writing forces act on atoms
        if task.arts.atomic_forces_status == True:
            task.arts.write_atomic_forces(fout)
        # =========================

    # gen pbs script
    with open("relax-cubic.pbs", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("#PBS -N %s\n" % args.jobname)
        fout.write("#PBS -l nodes=%d:ppn=%d\n" % (args.nodes, args.ppn))
        fout.write("\n")
        fout.write("cd $PBS_O_WORKDIR\n")
        fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")
        #fout.write("mpirun -np $NP -machinefile $PBS_NODEFILE %s < %s > %s\n" % (cmd, inpname, output))

        a = task.arts.xyz.cell[0][0]

        fout.write("for a in `seq -w %f %f %f`\n" % (a-args.na/2*args.stepa, args.stepa, a+args.na/2*args.stepa))
        fout.write("do\n")
        fout.write("  cp relax.in.template relax-${a}.in\n")
        fout.write("  cat >> relax-${a}.in <<EOF\n")
        fout.write("\n")
        fout.write("CELL_PARAMETERS angstrom\n")
        fout.write("${a} 0.000000 0.000000\n")
        fout.write("0.000000 ${a} 0.000000\n")
        fout.write("0.000000 0.000000 ${a}\n")
        fout.write("EOF\n")
        fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE pw.x < relax-${a}.in > relax-${a}.out\n")
        fout.write("done\n")

    # generate result analysis script
    os.system("mkdir -p post-processing")

    with open("post-processing/get_energy.sh", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("cat > energy-latconst.data <<EOF\n")
        fout.write("# format: a energy(Ry)\n")
        fout.write("EOF\n")
        fout.write("for a in `seq -w %f %f %f`\n" % (a-args.na/2*args.stepa, args.stepa, a+args.na/2*args.stepa))
        fout.write("do\n")
        fout.write("  energy=`cat ../relax-${a}.out | grep '!    total energy' | tail -1`\n")
        fout.write("  cat >> energy-latconst.data <<EOF\n")
        fout.write("${a} ${energy:32:-2}\n")
        fout.write("EOF\n")
        fout.write("done\n")
        fout.write("cat > energy-latconst.gp<<EOF\n")
        fout.write("set term gif\n")
        fout.write("set output 'energy-latconst.gif'\n")
        fout.write("set title Energy Latconst\n")
        fout.write("set xlabel 'latconst(a)'\n")
        fout.write("set ylabel 'Energy'\n")
        fout.write("plot 'energy-latconst.data' w l\n")
        fout.write("EOF\n")
    #os.system("cd post-processing; bash get_energy.sh; cd ../")
    os.chdir("../")

    server_handle(auto=args.auto, directory=args.directory, jobfilebase="relax-cubic", server=args.server)
