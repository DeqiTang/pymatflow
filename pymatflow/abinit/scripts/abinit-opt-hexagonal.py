#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import shutil
import argparse

from pymatflow.remote.server import server_handle
from pymatflow.abinit.opt import opt_run

"""
usage:
requirements:
    optcell should always be 0
    task.input.to_input() must make sure the acell is alway 1 1 1 angstrom
    cell is specified using acell and rprim
    and the cell parameter is right behind rprim in the next line
    like this:
    rprim
    3.192238 0.000000 0.000000
   -1.596119 2.764559 0.000000
    0.000000 0.000000 7.000000
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-abinit-opt-hexagonal",
            help="Directory to do the optimization calculation")

    parser.add_argument("-f", "--file", type=str,
            help="The xyz structure file with second line specifying cell parameters")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--ecut", type=int, default=15,
            help="Kinetic energy cutoff for wave functions in unit of Hartree, default value: 15 Hartree. for more information refer to https://docs.abinit.org/variables/basic/#ecut")

    parser.add_argument("--ixc", type=int, default=11,
            choices=[1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 40, 41, 42],
            help="type of exchage-correlation functional. fore more information, refer to https://docs.abinit.org/variables/basic/#ixc for more information")

    parser.add_argument("--kptopt", type=int, default=1,
            help="Kpoints Generation scheme option: 0, 1, 2, 3, 4 or a negative value. for more information, refer to https://docs.abinit.org/variables/basic/#kptopt")

    parser.add_argument("--ngkpt", nargs="+", type=int,
            default=[1, 1, 1],
            help="number of grid points for kpoints generation. for more information, refer to https://docs.abinit.org/variables/basic/#ngkpt")

    # vdw related parameters
    parser.add_argument("--vdw-xc", type=int, default=None,
            choices=[0, 1, 2, 5, 6, 7, 10, 11, 14],
            help="Van Der Waals exchange-correlation functional. 0: no correction, 1: vdW-DF1, 2: vdW-DF2, 5: DFT-D2, 6: DFT-D3, 7: DFT-D3(BJ). for more information, refer to https://docs.abinit.org/variables/vdw/#vdw_xc")

    parser.add_argument("--vdw-tol", type=float,
            default=None,
            help="Van Der Waals tolerance, only work when vdw_xc == 5 or 6 or 7. to be included in the potential a pair of atom must have contribution to the energy larger than vdw_tol. default value is 1.0e-10. fore more information, refer to https://docs.abinit.org/variables/vdw/#vdw_tol")

    # -----------------------------------------------------------
    #                        ions moving related parameters
    # -----------------------------------------------------------

    parser.add_argument("--ionmov", type=int, default=3,
            choices=[2, 3, 4, 5],
            help="type of ionmov algorithm. fore more information, refer to https://docs.abinit.org/variables/rlx/#ionmov")

    parser.add_argument("--optcell", type=int,
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            default=0,
            help="whether to optimize the cell shape and dimension. fore more information, refer to https://docs.abinit.org/variables/rlx/#optcell")

    parser.add_argument("--ecutsm", type=int, default=None,
            help="when optcell != 0, must specify encutsm larser than zero. for more information refer to https://docs.abinit.org/variables/rlx/#ecutsm")


    # ------------------------------------------------
    # na stepa nc stepc
    # ------------------------------------------------
    parser.add_argument("--na", type=int, default=10,
            help="number of a to run")
    parser.add_argument("--stepa", type=float, default=0.05,
            help="step of a in unit of Angstrom")
    parser.add_argument("--nc", type=int, default=10,
            help="number of c to run")
    parser.add_argument("--stepc", type=float, default=0.05,
            help="step of c in unit of Angstrom")

    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="opt-cubic",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================
    args = parser.parse_args()

    params = {}
    kpoints = {}

    params["ecut"] = args.ecut
    params["ixc"] = args.ixc
    params["vdw_xc"] = args.vdw_xc
    params["vdw_tol"] = args.vdw_tol

    params["kptopt"] = args.kptopt
    params["ngkpt"] = args.ngkpt

    params["ionmov"] = args.ionmov
    params["optcell"] = args.optcell
    params["ecutsm"] = args.ecutsm

    task = opt_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_kpoints(kpoints=kpoints)
    #task.optimize(directory=args.directory, mpi=args.mpi, runopt=args.runopt)
    #server_handle(auto=args.auto, directory=args.directory, jobfilebase="optimization", server=args.server)

    if os.path.exists(args.directory):
        shutil.rmtree(args.directory)
    os.mkdir(args.directory)
    os.system("cp *.psp8 %s/" % args.directory)
    os.system("cp *.GGA_PBE-JTH.xml %s/" % args.directory)
    os.system("cp %s %s/" % (task.input.system.xyz.file, args.directory))


    #task.poscar.to_poscar(os.path.join(args.directory, "POSCAR"))

    os.chdir(args.directory)
    # gen pbs script
    with open("opt-hexagonal.pbs", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("#PBS -N %s\n" % args.jobname)
        fout.write("#PBS -l nodes=%d:ppn=%d\n" % (args.nodes, args.ppn))
        fout.write("\n")
        fout.write("cd $PBS_O_WORKDIR\n")
        fout.write("cat > optimization.in<<EOF\n")
        #task.input.to_input(fout)
        fout.write(task.input.to_string())
        fout.write("EOF\n")
        fout.write("cat > optimization.files<<EOF\n")
        #task.files.name = "optimization.files"
        task.files.main_in = "optimization.in"
        task.files.main_out = "optimization.out"
        task.files.wavefunc_in = "optimization-i"
        task.files.wavefunc_out = "optimization-o"
        task.files.tmp = "tmp"
        #task.files.to_files(fout, task.input.system)
        fout.write(task.files.to_string(system=task.input.system))
        fout.write("EOF\n")
        fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

        a = task.input.system.xyz.cell[0][0]

        fout.write("v11=%f\n" % task.input.system.xyz.cell[0][0])
        fout.write("v12=%f\n" % task.input.system.xyz.cell[0][1])
        fout.write("v13=%f\n" % task.input.system.xyz.cell[0][2])
        fout.write("v21=%f\n" % task.input.system.xyz.cell[1][0])
        fout.write("v22=%f\n" % task.input.system.xyz.cell[1][1])
        fout.write("v23=%f\n" % task.input.system.xyz.cell[1][2])
        fout.write("v31=%f\n" % task.input.system.xyz.cell[2][0])
        fout.write("v32=%f\n" % task.input.system.xyz.cell[2][1])
        fout.write("v33=%f\n" % task.input.system.xyz.cell[2][2])

        fout.write("rprim_line=`cat optimization.in | grep -n \'rprim\' | cut -d \":\" -f 1`\n")
        fout.write("after_rprim_cell_line=`echo \"${rprim_line} + 4\" | bc`\n")
        if args.na >= 2:
            # a is optimized
            fout.write("for a in `seq -w %f %f %f`\n" % (a-args.na/2*args.stepa, args.stepa, a+args.na/2*args.stepa))
            fout.write("do\n")
            if args.nc >= 2:
                # optimize both a and c
                fout.write("for c in `seq -w %f %f %f`\n" % (c-args.nc/2*args.stepc, args.stepc, c+args.nc/2*args.stepc))
                fout.write("do\n")
                fout.write("  mkdir relax-${a}-${c}\n")
                fout.write("  cp  optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}-${c}/\n")
                fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}-${c}/optimization.in\n")
                fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                # here with the usage of length and scale in bs processing, we can make sure that number like '.123' will be correctly
                # set as '0.123', namely the ommited 0 by bs by default is not ommited now!
                fout.write("  cat >> relax-${a}-${c}/optimization.in<<EOF\n")
                fout.write("${a} 0.000000 0.000000\n")
                fout.write("${vec21} ${vec22} 0.000000\n")
                fout.write("0.000000 0.000000 ${c}\n")
                fout.write("EOF\n")
                fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}-${c}/optimization.in\n")
                fout.write("  cd relax-${a}-${c}/\n")
                fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE abinit < optimization.files\n")
                fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # only optimize a
                fout.write("  mkdir relax-${a}\n")
                fout.write("  cp optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${a}/\n")
                fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${a}/optimization.in\n")
                fout.write("  vec21=`echo \"scale=6; result=${v21} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                fout.write("  vec22=`echo \"scale=6; result=${v22} * ${a} / ${v11}; if (length(result)==scale(result)) print 0; print result\" | bc`\n")
                fout.write("  cat >> relax-${a}/optimization.in<<EOF\n")
                fout.write("${a} 0.000000 0.000000\n")
                fout.write("${vec21} ${vec22} 0.000000\n")
                fout.write("0.000000 0.000000 ${v33}\n")
                fout.write("EOF\n")
                fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${a}/optimization.in\n")
                fout.write("  cd relax-${a}/\n")
                fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE abinit < optimization.files\n")
                fout.write("  cd ../\n")
            fout.write("done\n")
        else:
            # a is not optimized
            if args.nc >= 2:
                # only optimize c
                fout.write("for c in `seq -w %f %f %f`\n" % (c-args.nc/2*args.stepc, args.stepc, c+args.nc/2*args.stepc))
                fout.write("do\n")
                fout.write("  mkdir relax-${c}\n")
                fout.write("  cp optimization.files *.psp8 *.GGA_PBE-JTH.xml relax-${c}/\n")
                fout.write("  cat optimization.in | head -n +${rprim_line} > relax-${c}/optimization.in\n")
                fout.write("  cat >> relax-${c}/optimization.in<<EOF\n")
                fout.write("${v11} 0.000000 0.000000\n")
                fout.write("${v21} ${v22} 0.000000\n")
                fout.write("0.000000 0.000000 ${c}\n")
                fout.write("EOF\n")
                fout.write("  cat optimization.in | tail -n +${after_rprim_cell_line} >> relax-${c}/optimization.in\n")
                fout.write("  cd relax-${c}/\n")
                fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE abinit < optimization.files\n")
                fout.write("  cd ../\n")
                fout.write("done\n")
            else:
                # neither a or c is optimized
                pass

    # generate result analysis script
    os.system("mkdir -p post-processing")

    with open("post-processing/get_energy.sh", 'w') as fout:
        fout.write("#!/bin/bash\n")
        # the comment
        if args.na >= 2 and args.nc >= 2:
            fout.write("cat > energy-latconst.data <<EOF\n")
            fout.write("# format: a c energy(eV)\n")
            fout.write("EOF\n")
        if args.na >= 2 and args.nc < 2:
            fout.write("cat > energy-latconst.data <<EOF\n")
            fout.write("# format: a energy(eV)\n")
            fout.write("EOF\n")
        if args.na < 2 and args.nc >= 2:
            fout.write("cat > energy-latconst.data <<EOF\n")
            fout.write("# format: c energy(eV)\n")
            fout.write("EOF\n")
        # end
        if args.na >= 2:
            fout.write("for a in `seq -w %f %f %f`\n" % (a-args.na/2*args.stepa, args.stepa, a+args.na/2*args.stepa))
            fout.write("do\n")
            if args.nc >= 2:
                fout.write("for c in `seq -w %f %f %f`\n" % (c-args.nc/2*args.stepc, args.stepc, c+args.nc/2*args.stepc))
                fout.write("do\n")
                fout.write("  energy=`cat ../relax-${a}-%{c}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % task.files.main_out)
                fout.write("  cat >> energy-latconst.data <<EOF\n")
                fout.write("${a} ${c} ${energy}\n")
                fout.write("EOF\n")
                fout.write("done\n")
                fout.write("done\n")

                fout.write("cat > energy-latconst.gp<<EOF\n")
                fout.write("set term gif\n")
                fout.write("set output 'energy-latconst.gif'\n")
                fout.write("set title Energy Latconst\n")
                fout.write("set xlabel 'latconst(a)'\n")
                fout.write("set ylabel 'latconst(c)'\n")
                fout.write("set zlabel 'Energy'\n")
                fout.write("splot 'energy-latconst.data'\n")
                fout.write("EOF\n")
            else:
                fout.write("  energy=`cat ../relax-${a}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % task.files.main_out)
                fout.write("  cat >> energy-latconst.data <<EOF\n")
                fout.write("${a} ${energy}\n")
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
        else:
            if args.nc >= 2:
                fout.write("for c in `seq -w %f %f %f`\n" % (c-args.nc/2*args.stepc, args.stepc, c+args.nc/2*args.stepc))
                fout.write("do\n")
                fout.write("  energy=`cat ../relax-${c}/%s | grep 'Etotal=' | tail -1 | cut -d \"=\" -f 1`\n" % task.files.main_out)
                fout.write("  cat >> energy-latconst.data <<EOF\n")
                fout.write("${c} ${energy}\n")
                fout.write("EOF\n")
                fout.write("done\n")
                fout.write("cat > energy-latconst.gp<<EOF\n")
                fout.write("set term gif\n")
                fout.write("set output 'energy-latconst.gif'\n")
                fout.write("set title Energy Latconst\n")
                fout.write("set xlabel 'latconst(c)'\n")
                fout.write("set ylabel 'Energy'\n")
                fout.write("plot 'energy-latconst.data' w l\n")
                fout.write("EOF\n")
            else:
                # nothing to do
                pass
    #os.system("cd post-processing; bash get_energy.sh; cd ../")
    os.chdir("../")


    server_handle(auto=args.auto, directory=args.directory, jobfilebase="opt-hexagonal", server=args.server)
