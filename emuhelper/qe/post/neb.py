#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys

class neb_post:
    def __init__(self):
        pass

    def min_energy_path_gp(self, directory="tmp-qe-neb", nebint='pwscf.int', nebdat='pwscf.dat', inpname="min-energy-path.gp", runopt="gen"):
        #
        # first check whether there is a previous neb running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("min_energy_path_gp plotting:\n")
            print("  directory of previous neb calculattion not found!\n")
            sys.exit(1)
        if runopt == "gen" or runopt == "genrun":
            with open(os.path.join(directory, inpname), 'w') as fout:
                fout.write("set term postscript enhanced\n")
                fout.write("set output 'min-energy-path.eps'\n")
                fout.write("set title 'Minmum Energy Path'\n")
                fout.write("set xlabel 'Reaction coordinate / arb. u.'\n")
                fout.write("set ylabel 'E - E_{IS} / eV'\n")
                fout.write("set format y '%.2f'\n")
                fout.write("set grid xtics ytics\n")
                fout.write("set xzeroaxis lt -1\n")
                fout.write("plot  [0:1][:] \\\n")
                fout.write("    '%s' notitle w l lt 2 lw 4, \\\n" % nebint)
                fout.write("    '%s' notitle w points lt 1 pt 7 ps 1.5\n" % nebdat)
                #fout.write("pause -1\n")

        if runopt == "run" or runopt == "genrun":
            os.chdir(directory)
            os.system("gnuplot %s" % inpname)
            os.chdir("../")


