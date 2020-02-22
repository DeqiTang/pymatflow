#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import datetime
import subprocess
import matplotlib.pyplot as plt


class md_post:
    """
    """
    def __init__(self, outputfile):
        """
        """
        self.job_completed = None # judge whether the calculation is finished
        self.md_params = {}
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

        self.get_md_params_and_run_info()

    def get_md_params_and_run_info(self):
        """
        """
        self.run_info["total-energies"] = []
        self.run_info["scf-iterations-converged"] = []
        self.md_params["VariableCell"] = "Unknown"
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
                self.md_params["VariableCell"] = line.split()[1]

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

    def markdown_report(self, md="MolecularDynamicsReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 分子动力学实验统计\n")
            fout.write("是否进行变胞分子动力学: %s\n" % self.md_params["VariableCell"])
            fout.write("分子动力学是否结束: %s\n" % str(self.job_completed))
            fout.write("## 分子动力学参数\n")
            for item in self.md_params:
                fout.write("- %s: %s\n" % (item, str(self.md_params[item])))
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
                fout.write("  - job is not finished yet, but it start at %s\n" % start)
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
        self.markdown_report("MolecularDynamicsReport.md")


