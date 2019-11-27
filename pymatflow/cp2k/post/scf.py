#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import datetime
import matplotlib.pyplot as plt


class scf_post:
    """
    """
    def __init__(self, output):
        """
        output is the output file of scf run
        """
        self.file = output
        self.scf_params = {}
        self.run_info = {}

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_info()

    def get_info(self):
        """
        get the general information of scf run from scf run output file
        which is now stored in self.lines
        """
        
        self.get_scf_params_and_run_info()
    # 
    def get_scf_params_and_run_info(self):
        """
        """
        #self.run_info["total-scf-energies"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "****" and line.split()[1] == "****" and line.split()[5] == "STARTED":
                self.run_info["start-time"] = line.split("\n")[0]
            if line.split()[0] == "****" and line.split()[1] == "****" and line.split()[5] == "ENDED":
                self.run_info["stop-time"] = line.split("\n")[0]
            if line.split()[0] == "eps_scf:":
                self.scf_params["EPS_SCF"] = float(line.split()[1])
            if line.split()[0] == "Mixing" and line.split()[1] == 'method:':
                self.scf_params["MIXING"] = line.split()[2]
            if line.split()[0] == "added" and line.split()[1] == "MOs":
                self.scf_params["ADDED_MOS"] = int(line.split()[2])
            if line.split()[0] == "***" and line.split()[1] == "SCF" and line.split()[3] == "converged":
                self.run_info["scf-steps"] = int(line.split()[5])
            if line.split()[0] == "ENERGY|" and line.split()[4] == "QS":
                self.run_info["total-energy"] = float(line.split()[8])
            if line.split()[0] ==  "Fermi" and line.split()[1] == "energy:":
                self.run_info["fermi-energy"] = float(line.split()[2])
    

    def plot_run_info(self):
        """
        """
        #plt.plot(self.run_info["total-scf-energies"])
        #plt.title("Total energy per scf step")
        #plt.xlabel("Scf step")
        #plt.ylabel("Total energy")
        #plt.tight_layout()
        #plt.savefig("total-energy-per-scf-step.png")
        #plt.close()
        pass


    def markdown_report(self, md="SCFReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# SCF实验统计\n")
            fout.write("## SCF参数\n")
            for item in self.scf_params:
                fout.write("- %s: %s\n" % (item, str(self.scf_params[item])))
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

    def export(self):
        self.plot_run_info()
        self.markdown_report("SCFReport.md")
