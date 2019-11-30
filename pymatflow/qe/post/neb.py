#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import datetime
import os
import sys

class neb_post:
    def __init__(self):
        self.neb_params = {}
        self.run_info = {}

    def get_info(self, nebout="neb.out"):
        with open(nebout, 'r') as fin:
            self.lines = fin.readlines()
            
        self.run_info["activation energy(->)"] = []
        self.run_info["activation eenrgy(<-)"] = []
        self.run_info["climbing image"] = []
        self.run_info["path length"] = []
        self.run_info["inter image distance"] = []
       
        for line in self.lines:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "Program" and line.split()[1] == "NEB" and line.split()[3] == "starts":
                self.run_info["start-time"] = line.split("\n")[0]
            if line.split()[0] == "This" and line.split()[1] == "run" and line.split()[3] == "terminated":
                self.run_info["stop-time"] = line.split("\n")[0]
            if line.split()[0] == "Program" and line.split()[3] == "starts":
                self.run_info["neb calculation starts at: "] = line.split()[5] + "-" + line.split()[7]
            if line.split()[0] == "This" and line.split()[3] == "terminated":
                self.run_info["neb calculation terminates at: "] = line.split()[6] + "-" + line.split()[5]
            if line.split()[0] == "Parallel"and line.split()[6] == "processors":
                self.run_info["processors"] = line.split()[5]
            if line.split()[0] == "MPI" and line.split()[5] == "nodes":
                self.run_info["nodes"] = line.split()[4]
            if line.split()[0] == "initial" and line.split()[1] == "path":
                self.run_info["initial-path-length(bohr)"] = line.split()[4]
            if line.split()[0] == "initial" and line.split()[1] == "inter-image":
                self.run_info["initial inter-image distance(bohr)"] = line.split()[4]
            if line.split()[0] == "string_method":
                self.neb_params["string_method"] = line.split()[2]
            if line.split()[0] == "opt_scheme":
                self.neb_params["opt_scheme"] = line.split()[2]
            if line.split()[0] == "num_of_images":
                self.neb_params["num_of_images"] = line.split()[2]
            if line.split()[0] == "nstep_path":
                self.neb_params["nstep_path"] = line.split()[2]
            if line.split()[0] == "CI_scheme":
                self.neb_params["CI_scheme"] = line.split()[2]
            if line.split()[0] == "first_last_opt":
                self.neb_params["first_last_opt"] = line.split()[2]
            if line.split()[0] == "use_freezing":
                self.neb_params["use_freezing"] = line.split()[2]
            if line.split()[0] == "ds":
                self.neb_params["ds(a.u.)"] = line.split()[2]
            if line.split()[0] == "k_max":
                self.neb_params["k_max(a.u.)"] = line.split()[2]
            if line.split()[0] == "k_min":
                self.neb_params["k_min(a.u.)"] = line.split()[2]
            if line.split()[0] == "suggested" and line.split()[1] == "k_max":
                self.neb_params["suggested k_max(a.u.)"] = line.split()[3]
            if line.split()[0] == "suggested" and line.split()[1] == "k_min":
                self.neb_params["suggested k_min(a.u.)"] = line.split()[3]
            if line.split()[0] == "path_thr":
                self.neb_params["path_thr(eV/A)"] = line.split()[2]
            if line.split()[0] == "activation" and line.split()[2] == "(->)":
                self.run_info["activation energy(->)"].append(float(line.split()[4]))
            if line.split()[0] == "activation" and line.split()[2] == "(<-)":
                self.run_info["activation eenrgy(<-)"].append(float(line.split()[4]))
            if line.split()[0] == "climbing" and line.split()[1] == "image":
                self.run_info["climbing image"].append(int(line.split()[3]))
            if line.split()[0] == "path" and line.split()[1] == "length":
                self.run_info["path length"].append(float(line.split()[3]))
            if line.split()[0] == "inter-image" and line.split()[1] == "distance":
                self.run_info["inter image distance"].append(float(line.split()[3]))



    def min_energy_path_gp(self, nebint='pwscf.int', nebdat='pwscf.dat', inpname="min-energy-path.gp", runopt="gen"):

        if runopt == "gen" or runopt == "genrun":
            with open(inpname, 'w') as fout:
                #fout.write("set term postscript enhanced\n")
                fout.write("set term gif\n")
                #fout.write("set output 'min-energy-path.eps'\n")
                fout.write("set output 'min-energy-path.gif'\n")
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
            os.system("gnuplot %s" % inpname)

    def md_report(self, md):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# Neb实验统计\n")
            fout.write("## neb参数\n")
            for item in self.neb_params:
                fout.write("- %s: %s\n" % (item, str(self.neb_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            # depending on the value of seconds in the time string, there are two situations:
            # when second is smaller than 10  it will be divided from xx:xx x, and when second
            # is larger or equal to 10 it will be together(xx:xx:xx).
            # so we have to preprocess it to build the right time string to pass into 
            # datetime.datetime.strptime()
            if len(self.run_info["start-time"].split()) == 8:
                start_str = self.run_info["start-time"].split()[5]+"-"+self.run_info["start-time"].split()[7]
            elif len(self.run_info["start-time"].split()) == 9:
                start_str = self.run_info["start-time"].split()[5]+"-"+self.run_info["start-time"].split()[7]+self.run_info["start-time"].split()[8]
            else:
                print("===============================================\n")
                print("                  Warning !!!\n")
                print("===============================================\n")
                print("qe.post.neb.markdown_report:\n")
                print("failed to parse start-time string\n")
                sys.exit(1)
            if len(self.run_info["stop-time"].split()) == 7:
                stop_str = self.run_info["stop-time"].split()[6]+"-"+self.run_info["stop-time"].split()[5]
            elif len(self.run_info["stop-time"].split()) == 8:
                stop_str = self.run_info["stop-time"].split()[7]+"-"+self.run_info["stop-time"].split()[5]+self.run_info["stop-time"].split()[6]
            else:
                print("===============================================\n")
                print("                  Warning !!!\n")
                print("===============================================\n")
                print("qe.post.neb.markdown_report:\n")
                print("failed to parse stop-time string\n")
                sys.exit(1)
            start = datetime.datetime.strptime(start_str, "%d%b%Y-%H:%M:%S")
            stop = datetime.datetime.strptime(stop_str, "%d%b%Y-%H:%M:%S")
            delta_t = stop -start
            fout.write("- Time consuming:\n")
            fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("Min energy path\n")
            fout.write("![min energy path](min-energy-path.gif)\n")
            fout.write("\n")

    def export(self, directory="tmp-qe-neb", nebint="pwscf.int", nebdat="pwscf.dat", nebout="neb.out", inpname="min-energy-path.gp", md="neb-report.md"):
        #
        # first check whether there is a previous neb running
        if not os.path.exists(directory):
            print("===================================================\n")
            print("                 Warning !!!\n")
            print("===================================================\n")
            print("min_energy_path_gp plotting:\n")
            print("  directory of previous neb calculattion not found!\n")
            sys.exit(1)       
        os.chdir(directory)
        self.get_info(nebout=nebout)
        self.min_energy_path_gp(nebint=nebint, nebdat=nebdat, inpname=inpname, runopt="genrun")
        self.md_report(md=md)
        os.chdir("../")
