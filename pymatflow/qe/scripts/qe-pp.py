#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse
from pymatflow.qe.static import static_run


"""
usage:
"""

inputpp = {}
plotpp = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-qe-static",
            help="directory of the static running")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz file name")

    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--auto", type=int, default=3,
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing in remote server, 3: pymatflow used in server with direct submit, in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    parser.add_argument("--plot-num", type=int, nargs="+", default=[0],
            choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18, 19, 20, 21],
            help="""
                type of analysis stored in the filplot file for later plot, 0: electron-pseudo-charge-density,
                    1: total-potential,
                    2: local-ionic-potential,
                    3: ldos,
                    4: local-density-of-electronic-entropy,
                    5: stm,
                    6: spin-polar,
                    7: molecular-orbitals,
                    8: electron-local-function,
                    9: charge-density-minus-superposition-of-atomic-densities,
                    10: ILDOS,
                    11: v_bare+v_H-potential,
                    12: sawtooth-electric-field-potential,
                    13: nocollinear-magnetization,
                    17: all-electron-charge-density-paw-only,
                    18: exchage-correlation-magnetic-field-noncollinear-case,
                    19: reduced-density-gradient,
                    20: product-of-charge-density-with-hessian,
                    21: all-electron-density-paw-only,""")

    parser.add_argument("--iflag", type=int,
            default=3,
            choices=[0, 1, 2, 3, 4],
            help="dimension of the plot. 0: 1D plot of the spherical average, 1: 1D plot, 2: 2D plot, 3: 3D plot, 4: 2D polar plot on a sphere")

    parser.add_argument("--output-format", type=int, default=5,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="output file format for visualization. 0: gnuplot(1D), 1: no longer supported, 2: plotrho(2D), 3: XCRYSDEN(2d), 4: no longer supported, 5: XCRYSDEN(3D), 6: gaussian cube(3D), 7: gnuplot(2D)")

    # -----------------------------------------------------------------
    #                       run params
    # -----------------------------------------------------------------

    parser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    parser.add_argument("--jobname", type=str, default="pp.x",
            help="jobname on the pbs server")

    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()
    xyzfile = args.file

    inputpp["plot_num"] = args.plot_num
    plotpp["iflag"] = args.iflag
    plotpp["output_format"] = args.output_format

    task = static_run()
    task.get_xyz(xyzfile)
    task.set_pp(inputpp=inputpp, plotpp=plotpp)
    task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
    task.pp(directory=args.directory, runopt=args.runopt, auto=args.auto)
