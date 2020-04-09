#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import datetime
import subprocess
import matplotlib.pyplot as plt



class opt_out:
    """
    """
    def __init__(self):
        """
        """
        self.job_completed = None # judge whether the calculation is finished
        self.opt_params = {}
        self.run_info = {}



    def get_info(self, filepath):
        """
        filepath: output file of the optimization running
        """
        self.file = filepath
        with open(self.file, 'r') as fin:
            self.lines = fin.readlines()

        # check whether calculation is finished
        if len(self.lines[-1].split()) == 2 and self.lines[-1].split()[0] == "Job" and self.lines[-1].split()[1] == "completed":
            self.job_completed = True
        else:
            self.job_completed = False
        #
        # check whether successfully relaxed
        self.relaxed = False
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "outcoor:" and self.lines[i].split()[1] == "Relaxed":
                self.relaxed = True
                break
        #
        self.get_opt_params_and_run_info()

    def get_opt_params_and_run_info(self):
        """
        """
        self.run_info["total_energy_each_ion_step"] = []
        self.run_info["scf_iterations_each_ion_step_converged"] = []
        for line in self.lines:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == ">>" and line.split()[1] == "Start":
                self.run_info["start_time"] = line.split("\n")[0]
            if line.split()[0] == ">>" and line.split()[1] == "End":
                self.run_info["stop_time"] = line.split("\n")[0]
            if line.split()[0] == "siesta:":
                if len(line.split()) ==  4 and line.split()[1] == "E_KS(eV)":
                    self.run_info["total_energy_each_ion_step"].append(float(line.split()[3]))
            if line.split()[0] == "SCF" and line.split()[1] == "cycle" and line.split()[2] == "converged":
                self.run_info["scf_iterations_each_ion_step_converged"].append(int(line.split()[4]))
            if line.split()[0] == "MD.VariableCell":
                self.opt_params["VariableCell"] = line.split()[1]

    def plot_info(self):
        #plot_run_info
        plt.plot(self.run_info["total_energy_each_ion_step"])
        plt.title("Total energy each ion step")
        plt.xlabel("Scf cycle")
        plt.ylabel("Total Energy (eV)")
        plt.tight_layout()
        plt.savefig("total-energy-each-ion-step.png")
        plt.close()

        plt.plot(self.run_info["scf_iterations_each_ion_step_converged"])
        plt.title("Iterations each ion step(converged)")
        plt.xlabel("Scf cycle(converged)")
        plt.ylabel("Scf iterations")
        plt.tight_layout()
        plt.savefig("iterations-each-ion-step-converged.png")
        plt.close()


    def markdown_report(self, md="opt-info.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("是否进行变胞优化: %s\n" % self.opt_params["VariableCell"])
            fout.write("几何优化任务是否结束:%s\n" % str(self.job_completed))
            if self.job_completed == True:
                fout.write("是否成功优化: %s\n" % str(self.relaxed))
            else:
                fout.write("是否成功优化: %s\n" % "运行未结束, 结果未知")
            fout.write("## 优化参数\n")
            for item in self.opt_params:
                fout.write("- %s: %s\n" % (item, str(self.opt_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            start = datetime.datetime.strptime(self.run_info["start_time"].split()[4]+"-"+self.run_info["start_time"].split()[5], "%d-%b-%Y-%H:%M:%S")
            if self.job_completed == True:
                stop = datetime.datetime.strptime(self.run_info["stop_time"].split()[4]+"-"+self.run_info["stop_time"].split()[5], "%d-%b-%Y-%H:%M:%S")
                delta_t = stop -start
            fout.write("- Time consuming:\n")
            if self.job_completed == True:
                fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            else:
                fout.write("  - job is not finished yet, but it starts at %s\n" % start)
            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("Iterations per SCF\n")
            fout.write("![Iterations per SCF](iterations-each-ion-step-converged.png)\n")

            fout.write("Total energies per SCF\n")
            fout.write("![Total energies per SCF](total-energy-each-ion-step.png)\n")

    def export(self, directory):
        """
        """
        os.chdir(directory)
        os.system("mkdir -p post-processing")
        os.chdir("post-processing")
        self.plot_info()
        self.markdown_report()
        os.chdir("../../")



class opt_post:
    """
    """
    def __init__(self, outputfile):
        """
        """
        self.job_completed = None # judge whether the calculation is finished
        self.opt_params = {}
        self.run_info = {}

        with open(outputfile, 'r') as fin:
            self.lines = fin.readlines()

        self.get_info()

    def get_info(self):
        """
        """
        # check whether calculation is finished
        if len(self.lines[-1].split()) == 2 and self.lines[-1].split()[0] == "Job" and self.lines[-1].split()[1] == "completed":
            self.job_completed = True
        else:
            self.job_completed = False
        #
        # check whether successfully relaxed
        self.relaxed = False
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "outcoor:" and self.lines[i].split()[1] == "Relaxed":
                self.relaxed = True
                break
        #
        self.get_opt_params_and_run_info()

    def get_opt_params_and_run_info(self):
        """
        """
        self.run_info["total-energies"] = []
        self.run_info["scf-iterations-converged"] = []
        for line in self.lines:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == ">>" and line.split()[1] == "Start":
                self.run_info["start-time"] = line.split("\n")[0]
            if line.split()[0] == ">>" and line.split()[1] == "End":
                self.run_info["stop-time"] = line.split("\n")[0]
            if line.split()[0] == "siesta:":
                if len(line.split()) ==  4 and line.split()[1] == "E_KS(eV)":
                    self.run_info["total-energies"].append(float(line.split()[3]))
            if line.split()[0] == "SCF" and line.split()[1] == "cycle" and line.split()[2] == "converged":
                self.run_info["scf-iterations-converged"].append(int(line.split()[4]))
            if line.split()[0] == "MD.VariableCell":
                self.opt_params["VariableCell"] = line.split()[1]

    def view_trajectory(self, trajfile="siesta.ANI"):
        subprocess.call(["xcrysden", "--xyz", trajfile])

    def plot_run_info(self):
        plt.plot(self.run_info["total-energies"])
        plt.title("Total energies per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Total Energies (eV)")
        plt.tight_layout()
        plt.savefig("total-energies-per-scf.png")
        plt.close()

        plt.plot(self.run_info["scf-iterations-converged"])
        plt.title("Iterations per SCF(converged)")
        plt.xlabel("Scf cycle(converged)")
        plt.ylabel("Scf iterations")
        plt.tight_layout()
        plt.savefig("iterations-per-scf-converged.png")
        plt.close()

    def markdown_report(self, md="OptimizationReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("是否进行变胞优化: %s\n" % self.opt_params["VariableCell"])
            fout.write("几何优化任务是否结束:%s\n" % str(self.job_completed))
            if self.job_completed == True:
                fout.write("是否成功优化: %s\n" % str(self.relaxed))
            else:
                fout.write("是否成功优化: %s\n" % "运行未结束, 结果未知")
            fout.write("## 优化参数\n")
            for item in self.opt_params:
                fout.write("- %s: %s\n" % (item, str(self.opt_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            start = datetime.datetime.strptime(self.run_info["start-time"].split()[4]+"-"+self.run_info["start-time"].split()[5], "%d-%b-%Y-%H:%M:%S")
            if self.job_completed == True:
                stop = datetime.datetime.strptime(self.run_info["stop-time"].split()[4]+"-"+self.run_info["stop-time"].split()[5], "%d-%b-%Y-%H:%M:%S")
                delta_t = stop -start
            fout.write("- Time consuming:\n")
            if self.job_completed == True:
                fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            else:
                fout.write("  - job is not finished yet, but it starts at %s\n" % start)
            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("Iterations per SCF\n")
            fout.write("![Iterations per SCF](iterations-per-scf-converged.png)\n")

            fout.write("Total energies per SCF\n")
            fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")

    def export(self):
        self.plot_run_info()
        self.markdown_report("OptimizationReport.md")
