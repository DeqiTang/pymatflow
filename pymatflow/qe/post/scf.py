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
        self.run_info["total-scf-energies"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "Program" and line.split()[1] == "PWSCF" and line.split()[3] == "starts":
                self.run_info["start-time"] = line.split("\n")[0]
            if line.split()[0] == "This" and line.split()[1] == "run" and line.split()[3] == "terminated":
                self.run_info["stop-time"] = line.split("\n")[0]
            if line.split()[0] == "kinetic-energy":
                self.scf_params["ecutwfc"] = int(float(line.split()[3]))
            if line.split()[0] == "convergence threshold":
                self.scf_params["conv_thr"] = float(line.split()[3])
            if line.split()[0] == "mixing" and line.split()[1] == 'beta':
                self.scf_params["mixing_beta"] = float(line.split()[3])
            if line.split()[0] == "number" and line.split()[2] == 'k':
                if line.split()[5] == "(tetrahedron":
                    self.scf_params["degauss"] = "tetrahedron method: degauss not needed"
                else:
                    self.scf_params["degauss"] = float(line.split()[9])
                self.run_info["number-of-k-points"] = int(line.split()[4])
            if line.split()[0] == "convergence" and line.split()[3] == "achieved":
                self.run_info["iterations"] = int(line.split()[5])
            if line.split()[0] == "total" and line.split()[1] == "energy":
                self.run_info["total-scf-energies"].append(float(line.split()[3]))
            if line.split()[0] == "!" and line.split()[5] == "Ry":
                self.run_info["total-energy"] = float(line.split()[4])
            if line.split()[0] ==  "the" and line.split()[1] == "Fermi":
                self.run_info["fermi-energy"] = float(line.split()[4])
            if line.split()[0] == "Total" and line.split()[1] == "force":
                self.run_info["total-force"] = float(line.split()[3])
            if line.split()[0] == "number" and line.split()[2] == "electrons":
                self.scf_params["number-of-electrons"] = int(float(line.split()[4]))
    
        # note: at present len(self.run_info["total-scf-energies"]) = len(self.run_info["iterations"]) - 1
        # because the total energy of the last step is not printed in format like the previous scf step,
        # and it is printed as the '!    total energy              = ' where there is a "!" in the beginning
        # now we append the final scf step energy to self.run_info["total-scf-energies"]
        self.run_info["total-scf-energies"].append(self.run_info["total-energy"])

    def plot_run_info(self):
        """
        """
        plt.plot(self.run_info["total-scf-energies"])
        plt.title("Total energy per scf step")
        plt.xlabel("Scf step")
        plt.ylabel("Total energy")
        plt.tight_layout()
        plt.savefig("total-energy-per-scf-step.png")
        plt.close()



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
                print("qe.post.scf.markdown_report:\n")
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
                print("qe.post.scf.markdown_report:\n")
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
            fout.write("Total energy per scf step\n")
            fout.write("![Total energy per scf step](total-energy-per-scf-step.png)\n")

    def export(self):
        self.plot_run_info()
        self.markdown_report("SCFReport.md")