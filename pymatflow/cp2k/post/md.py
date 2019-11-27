#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import datetime
import matplotlib.pyplot as plt

class md_post:
    """
    """
    def __init__(self, output, run_type="MD"):
        self.file = output
        self.run_type = run_type
        self.md_params = {}
        self.run_info = {}

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_info()

    def get_info(self, output="md.out"):
        """
        get the general information of opt run from opt run output file
        which is now stored in self.lines
        """
        self.get_final_structure()

        self.get_md_params_and_run_info()

    def get_final_structure(self):
        pass

    def get_md_params_and_run_info(self):
        """
        """
        self.run_info["scf-steps-converged"] = []
        self.run_info["total-energies"] = []
        self.run_info["fermi-energies"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "****" and line.split()[1] == "****" and line.split()[5] == "STARTED":
                self.run_info["start-time"] = line.split("\n")[0]
            if line.split()[0] == "****" and line.split()[1] == "****" and line.split()[5] == "ENDED":
                self.run_info["stop-time"] = line.split("\n")[0]
            if line.split()[0] == "eps_scf:":
                self.md_params["EPS_SCF"] = float(line.split()[1])
            if line.split()[0] == "Mixing" and line.split()[1] == 'method:':
                self.md_params["MIXING"] = line.split()[2]
            if line.split()[0] == "added" and line.split()[1] == "MOs":
                self.md_params["ADDED_MOS"] = int(line.split()[2])
            # there might be scf cycles that are not converged
            # so the total scf cycles might be more than the len(self.run_info["scf-steps-converged"])
            if line.split()[0] == "***" and line.split()[1] == "SCF" and line.split()[3] == "converged":
                self.run_info["scf-steps-converged"].append(int(line.split()[5]))
            # in every scf cycle, whether it converged or not, the ENERGY| Total FORCE_EVAL ( QS  ) energy (a.u.):
            # will always be print out
            # so len(self.run_info["total-energies"]) equal to the actually times for scf calculation
            if line.split()[0] == "ENERGY|" and line.split()[4] == "QS":
                self.run_info["total-energies"].append(float(line.split()[8]))

        #self.run_info["scf-cycles"] = len(self.run_info["scf-steps"])
        # the total scf cycles might be more than the len(self.run_info["scf-steps"])
        # so we should use self.run_info["total-energies"] to get the total number of scf cycles
        self.run_info["scf-cycles-total"] = len(self.run_info["total-energies"])

    def plot_run_info(self):
        """
        """
        # this only plot the converged scf cycles.
        plt.plot(self.run_info["scf-steps-converged"])
        plt.title("Iterations per SCF(converged)")
        plt.xlabel("Scf cycle(converged)")
        plt.ylabel("Scf iterations")
        plt.tight_layout()
        plt.savefig("iterations-per-scf-converged.png")
        plt.close()

        plt.plot(self.run_info["total-energies"])
        plt.title("Total energies per SCF")
        plt.xlabel("Scf cycle")
        plt.ylabel("Total Energies (a.u.)")
        plt.tight_layout()
        plt.savefig("total-energies-per-scf.png")
        plt.close()


    def markdown_report(self, md="MolecularDynamicsReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 分子动力学实验统计\n")
            fout.write("分子动力学类型: %s\n" % self.run_type)
            fout.write("## 分子动力学参数\n")
            for item in self.md_params:
                fout.write("- %s: %s\n" % (item, str(self.md_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            # note that the accuracy for the seconds is not fully guranteed 
            # e.g. 2019-11-26 12:09:36.487 is read as 2019-11-26 12:09:36
            start = datetime.datetime.strptime(self.run_info["start-time"].split()[7]+"-"+self.run_info["start-time"].split()[8].split(".")[0], "%Y-%m-%d-%H:%M:%S")
            stop = datetime.datetime.strptime(self.run_info["stop-time"].split()[7]+"-"+self.run_info["stop-time"].split()[8].split(".")[0], "%Y-%m-%d-%H:%M:%S")
            delta_t = stop -start
            fout.write("- Time consuming:\n")
            fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("Iterations per SCF(converged)\n")
            fout.write("![Iterations per SCF](iterations-per-scf-converged.png)\n")
            
            fout.write("Total energies per SCF\n")
            fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")


    def export(self):
        self.plot_run_info()
        self.markdown_report("MolecularDynamicsReport.md")

