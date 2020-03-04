#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import shutil
import argparse

from pymatflow.remote.server import server_handle

from pymatflow.vasp.static import static_run

"""
usage:
"""

params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-vasp-efficiency-test",
            help="directory of the static running")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz file name")

    parser.add_argument("--runopt", type=str, default="gen",
            help="gen, run, or genrun")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--restart", type=int, default=0,
            choices=[0, 1],
            help="restart option. 0: new caculation; 1: restart caclculation")

    # --------------------------------------------------------
    #                   INCAR PARAMETERS
    # --------------------------------------------------------
    parser.add_argument("--prec", type=str, default="Normal",
            choices=["Normal", "Accurate", "A", "N"],
            help="PREC, default value: Normal")

    parser.add_argument("--encut", type=int, default=300,
            help="ENCUT, default value: 300 eV")

    parser.add_argument("--ediff", type=float, default=1.0e-4,
            help="EDIFF, default value: 1.0e-4")

    #parser.add_argument("--kpoints-mp", type=int, nargs="+",
    #        default=[1, 1, 1, 0, 0, 0],
    #        help="set kpoints like -k 1 1 1 0 0 0")
    parser.add_argument("--krange", type=int, nargs="+",
            default=[1, 4, 1],
            help="test k range")

    parser.add_argument("--ismear", type=int, default=0,
            help="smearing type(methfessel-paxton(>0), gaussian(0), fermi-dirac(-1), tetra(-4), tetra-bloch-dorrected(-5)), default: 0")

    parser.add_argument("--sigma", type=float, default=0.01,
            help="determines the width of the smearing in eV.")

    parser.add_argument("--ivdw", type=int, default=None,
            choices=[0, 11, 12, 21, 202, 4],
            help="IVDW = 0(no correction), 1(dft-d2), 11(dft-d3 Grimme), 12(dft-d3 Becke-Jonson), default: None which means 0, no correction")


    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh", "lsf_sz"],
            help="type of remote server, can be pbs or yh or lsf_sz")

    parser.add_argument("--jobname", type=str, default="efficiency-test",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to static_run setting
    # ==========================================================
    args = parser.parse_args()

    params["ENCUT"] = args.encut
    params["EDIFF"] = args.ediff
    params["ISMEAR"] = args.ismear
    params["SIGMA"] = args.sigma
    params["IVDW"] = args.ivdw



    task = static_run()
    task.get_xyz(args.file)
    task.set_params(params)
    #task.set_kpoints(kpoints_mp=args.kpoints_mp)
    #task.scf(directory=args.directory, runopt=args.runopt, restart=args.restart, mpi=args.mpi, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)

    if os.path.exists(args.directory):
        shutil.rmtree(args.directory)
    os.mkdir(args.directory)
    shutil.copyfile("POTCAR", os.path.join(args.directory, "POTCAR"))
    os.system("cp %s %s/" % (task.poscar.xyz.file, args.directory))

    #with open(os.path.join(args.directory, "INCAR"), 'w') as fout:
    #    task.incar.to_incar(fout)
    #task.kpoints.to_kpoints(os.path.join(args.directory, "KPOINTS"))
    task.poscar.to_poscar(os.path.join(args.directory, "POSCAR"))

    os.chdir(args.directory)

    # gen pbs script
    with open("efficiency-test-k.pbs", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("#PBS -N %s\n" % args.jobname)
        fout.write("#PBS -l nodes=%d:ppn=%d\n" % (args.nodes, args.ppn))
        fout.write("\n")
        fout.write("cd $PBS_O_WORKDIR\n")
        fout.write("cat > INCAR<<EOF\n")
        task.incar.to_incar(fout)
        fout.write("EOF\n")
        fout.write("NP=`cat $PBS_NODEFILE | wc -l`\n")

        fout.write("for k in `seq -w %d %d %d`\n" % (args.krange[0], args.krange[2], args.krange[1]))
        fout.write("do\n")
        fout.write("  mkdir -p kpoint-${k}\n")
        fout.write("  cp POTCAR INCAR POSCAR kpoint-${k}/\n")
        fout.write("  cat > kpoint-${k}/KPOINTS<<EOF\n")
        fout.write("K-POITNS automatic\n")
        fout.write("0\n")
        fout.write("Gamma\n")
        fout.write("${k} ${k} ${k}\n")
        fout.write("0 0 0\n")
        fout.write("EOF\n")
        fout.write("  cd kpoint-${k}\n")
        fout.write("  mpirun -np $NP -machinefile $PBS_NODEFILE -genv I_MPMI_FABRICS shm:tmi vasp_std\n")
        fout.write("  cd ../\n")
        fout.write("done\n")

    # result analysis script
    with open("get_info.sh", 'w') as fout:
        fout.write("#!/bin/bash\n")
        fout.write("cat > time.data<<EOF\n")
        fout.write("# format: kxkxk k cpu_time user_time system_time elapsed_time\n")
        fout.write("EOF\n")
        fout.write("for k in `seq -w %d %d %d`\n" % (args.krange[0], args.krange[2], args.krange[1]))
        fout.write("do\n")
        fout.write("  cpu_time=`grep time kpoint-${k}/OUTCAR | tail -n -4 | head -n 1`\n")
        fout.write("  user_time=`grep time kpoint-${k}/OUTCAR | tail -n -4 | head -n 2 | tail -n -1`\n")
        fout.write("  system_time=`grep time kpoint-${k}/OUTCAR | tail -n -4 | head -n 3 | tail -n -1`\n")
        fout.write("  elapsed_time=`grep time kpoint-${k}/OUTCAR | tail -n -4 | head -n 4 | tail -n -1`\n")
        fout.write("  cat >> time.data<<EOF\n")
        fout.write("${k}x${k}x${k} ${k} ${cpu_time:44} ${user_time:44} ${system_time:44} ${elapsed_time:44}\n")
        fout.write("EOF\n")
        fout.write("done\n")
        fout.write("\n")
        fout.write("cat > plot.gnuplot<<EOF\n")
        fout.write("set term png\n")
        fout.write("set output 'efficiency-k.png'\n")
        fout.write("set title 'Vasp Efficiency Test on Kpoints'\n")
        fout.write("set xlabel 'kpoints'\n")
        fout.write("set ylabel 'Time (sec)'\n")
        fout.write("set xtics(")
        for k in range(args.krange[0], args.krange[1], args.krange[2]):
            fout.write("'%dx%dx%d' %d, " % (k, k, k, k))
        fout.write("'%dx%dx%d' %d)\n" % (args.krange[1], args.krange[1], args.krange[1], args.krange[1]))
        fout.write("plot 'time.data' u 2:3 w l title 'CPU time',\\\\\n")
        fout.write("  'time.data' u 2:4 w l title 'USER time',\\\\\n")
        fout.write("  'time.data' u 2:5 w l title 'SYSTEM time',\\\\\n")
        fout.write("  'time.data' u 2:6 w l title 'ELAPSED time',\\\\\n")
        fout.write("EOF\n")
        fout.write("\n")
        fout.write("gnuplot plot.gnuplot\n")

    os.chdir("../")

    server_handle(auto=args.auto, directory=args.directory, jobfilebase="efficient-test-k", server=args.server)
