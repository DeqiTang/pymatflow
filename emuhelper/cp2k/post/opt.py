#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import matplotlib.pyplot as plt

class opt_post:
    """
    """
    def __init__(self, output, run_type):
        self.file = output
        self.run_type = run_type
        self.cell = None # optimized cell
        self.atoms = None # optimized atoms
        self.opt_params = {}
        self.run_info = {}

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_info()

    def get_info(self, output="geo-opt.out"):
        """
        get the general information of opt run from opt run output file
        which is now stored in self.lines
        """
        self.get_final_structure()

        self.get_opt_params_and_run_info()

    def get_final_structure(self):
        pass

    def get_opt_params_and_run_info(self):
        """
        """
        self.run_info["scf-steps"] = []
        self.run_info["total-energies"] = []
        self.run_info["fermi-energies"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "eps_scf:":
                self.opt_params["EPS_SCF"] = float(line.split()[1])
            if line.split()[0] == "Mixing" and line.split()[1] == 'method:':
                self.opt_params["MIXING"] = line.split()[2]
            if line.split()[0] == "added" and line.split()[1] == "MOs":
                self.opt_params["ADDED_MOS"] = int(line.split()[2])
            if line.split()[0] == "***" and line.split()[1] == "SCF" and line.split()[3] == "converged":
                self.run_info["scf-steps"].append(int(line.split()[5]))
            if line.split()[0] == "ENERGY|" and line.split()[4] == "QS":
                self.run_info["total-energies"].append(float(line.split()[8]))
            if line.split()[0] ==  "Fermi" and line.split()[1] == "energy:":
                self.run_info["fermi-energies"].append(float(line.split()[2]))

        self.run_info["scf-cycles"] = len(self.run_info["scf-steps"])

    def plot_run_info(self):
        """
        """
        plt.plot(self.run_info["scf-steps"])
        plt.title("Iterations per SCF")
        plt.xlabel("Scf cycle")
        plt.ylabel("Scf iterations")
        plt.tight_layout()
        plt.savefig("iterations-per-scf.png")
        plt.close()

        plt.plot(self.run_info["total-energies"])
        plt.title("Total energies per SCF")
        plt.xlabel("Scf cycle")
        plt.ylabel("Total Energies (a.u.)")
        plt.tight_layout()
        plt.savefig("total-energies-per-scf.png")
        plt.close()

        plt.plot(self.run_info["fermi-energies"])
        plt.title("Fermi energies per SCF")
        plt.xlabel("Scf cycle")
        plt.ylabel("Fermi energies ()")
        plt.tight_layout()
        plt.savefig("fermi-energies-per-scf.png")
        plt.close()

    def markdown_report(self, md="OptimizationReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("几何优化类型: %s\n" % self.run_type)
            fout.write("## 优化参数\n")
            for item in self.opt_params:
                fout.write("- %s: %s\n" % (item, str(self.opt_params[item])))
            fout.write("## 运行信息\n")
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("Iterations per SCF\n")
            fout.write("![Iterations per SCF](iterations-per-scf.png)\n")
            
            fout.write("Total energies per SCF\n")
            fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")

            fout.write("Fermi energies per SCF\n")
            fout.write("![Fermi energies per SCF](fermi-energies-per-scf.png)\n")


    def export(self):
        self.plot_run_info()
        self.markdown_report("OptimizationReport.md")
