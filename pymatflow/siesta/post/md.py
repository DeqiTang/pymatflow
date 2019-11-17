#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import matplotlib.pyplot as plt


class md_post:
    """
    """
    def __init__(self, outputfile):

        self.md_params = {}
        self.run_info = {}

        with open(outputfile, 'r') as fin:
            self.lines = fin.readlines()
    
    def process(self):
        self.get_md_params_and_run_info()
        self.plot_run_info()

    def get_md_params_and_run_info(self):
        """
        """
        self.run_info["total-energies"] = []
        self.run_info["scf-iterations-converged"] = []
        self.md_params["VariableCell"] = "Unknown"
        for line in self.lines:
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "siesta:":
                if len(line.split()) ==  4 and line.split()[1] == "E_KS(eV)":
                    self.run_info["total-energies"].append(float(line.split()[3]))
            if line.split()[0] == "SCF" and line.split()[1] == "cycle" and line.split()[2] == "converged":
                self.run_info["scf-iterations-converged"].append(int(line.split()[4]))
            if line.split()[0] == "MD.VariableCell":
                self.md_params["VariableCell"] = line.split()[1]


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
            fout.write("## 分子动力学参数\n")
            for item in self.md_params:
                fout.write("- %s: %s\n" % (item, str(self.md_params[item])))
            fout.write("## 运行信息\n")
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("Iterations per SCF\n")
            fout.write("![Iterations per SCF](iterations-per-scf-converged.png)\n")
            
            fout.write("Total energies per SCF\n")
            fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")

    def export(self):
        self.process()
        self.markdown_report("MolecularDynamicsReport.md")


