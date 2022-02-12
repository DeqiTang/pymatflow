# ==============================================================================
import os
import copy
import sys
import json
import datetime
import matplotlib.pyplot as plt

from pymatflow.base.atom import Atom
from pymatflow.base.xyz import BaseXyz

class ScfOut:
    """
    """
    def __init__(self):
        """
        output is the output file of scf run
        """
        self.file = None
        self.scf_params = {}
        self.run_info = {}

    def get_info(self, file):
        """
        get the general information of scf run from scf run output file
        which is now stored in self.lines
        """
        self.clean()

        self.file = file
        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_scf_params_and_run_info()

    def clean(self):
        self.file = None
        self.scf_params = {}
        self.run_info = {}

    #
    def get_scf_params_and_run_info(self):
        """
        self.run_info[]
            start_time: the task start time
            stop_time: the task stop time
            scf_energies: all the energies during the scf procedure
            #fermi_energy: fermi energy of the system (if output)

        """
        self.run_info["scf_energies"] = []

        for i in range(len(self.lines)):
            # if it is an empty line continue to next line
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == "****" and self.lines[i].split()[1] == "****" and self.lines[i].split()[5] == "STARTED":
                self.run_info["start_time"] = self.lines[i].split("\n")[0]
            elif self.lines[i].split()[0] == "****" and self.lines[i].split()[1] == "****" and self.lines[i].split()[5] == "ENDED":
                self.run_info["stop_time"] = self.lines[i].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Basis":
                self.scf_params["BASIS_SET_FILE_NAME"] = self.lines[i].split()[-1].split("\n")[0]
                self.scf_params["POTENTIAL_FILE_NAME"] = self.lines[i+1].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Project":
                self.scf_params["PROJECT_NAME"] = self.lines[i].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Run":
                self.scf_params["RUN_TYPE"] = self.lines[i].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Global":
                self.scf_params["PRINT_LEVEL"] = self.lines[i].split()[-1].split("\n")[0]
            elif self.lines[i].split()[0] == "GLOBAL|" and self.lines[i].split()[1] == "Total" and self.lines[i].split()[2] == "number":
                self.run_info["mpi_processes"] = int(self.lines[i].split()[-1].split("\n")[0])
                self.run_info["cpu_model"] = self.lines[i+3].split("name")[1].split("\n")[0]
            elif self.lines[i].split()[0] == "-" and self.lines[i].split()[1] == "Atoms:":
                self.run_info["n_atom"] = int(self.lines[i].split()[-1].split("\n")[0])
                self.run_info["n_shell"] = int(self.lines[i+2].split()[-1].split("\n")[0])
            elif self.lines[i].split()[0] == "max_scf:":
                self.scf_params["MAX_SCF"] = int(self.lines[i].split()[-1].split("\n")[0])
                self.scf_params["MAX_SCF_HISTORY"] = int(self.lines[i+1].split()[-1].split("\n")[0])
                self.scf_params["MAX_DIIS"] = int(self.lines[i+2].split()[-1].split("\n")[0])
            elif self.lines[i].split()[0] == "eps_scf:":
                self.scf_params["EPS_SCF"] = float(self.lines[i].split()[1])
                self.scf_params["EPS_SCF_HISTORY"] = float(self.lines[i+1].split()[1])
                self.scf_params["EPS_DIIS"] = float(self.lines[i+2].split()[1])
            elif self.lines[i].split()[0] == "Mixing" and self.lines[i].split()[1] == 'method:':
                self.scf_params["MIXING"] = self.lines[i].split()[2]
            elif self.lines[i].split()[0] == "added" and self.lines[i].split()[1] == "MOs":
                self.scf_params["ADDED_MOS"] = int(self.lines[i].split()[2])
            elif self.lines[i].split()[0] == "Number" and self.lines[i].split()[1] == "of" and self.lines[i].split()[2] == "electrons:":
                self.run_info["n_electrons"] = int(self.lines[i].split()[3])
                self.run_info["n_occcupied_orbital"] = int(self.lines[i+1].split()[4])
                self.run_info["n_molecular_orbital"] = int(self.lines[i+2].split()[4])
                self.run_info["n_orbital_function"] = int(self.lines[i+4].split()[4])
            elif self.lines[i].split()[0] == "Step" and self.lines[i].split()[1] == "Update" and self.lines[i].split()[-1] == "Change":
                self.run_info["scf_head_line"] = i
            elif self.lines[i].split()[0] == "***" and self.lines[i].split()[1] == "SCF" and self.lines[i].split()[2] == "run" and self.lines[i].split()[3] == "converged":
                self.run_info["scf_steps"] = int(self.lines[i].split()[5])
            elif self.lines[i].split()[0] ==  "Fermi" and self.lines[i].split()[1] == "energy:":
                self.run_info["fermi_energy"] = float(self.lines[i].split()[2])
            elif self.lines[i].split()[0] == "Total" and self.lines[i].split()[1] == "energy:":
                self.run_info["final_scf_energy"] = float(self.lines[i].split()[2])
            elif self.lines[i].split()[0] == "ENERGY|" and self.lines[i].split()[1] == "Total" and self.lines[i].split()[2] == "FORCE_EVAL":
                self.run_info["free_energy"] = float(self.lines[i].split()[8])
            else:
                pass
        # now we obtain the total energy for each dft scf step
        # note: cp2k will output the scf DFT energy for each scf step.
        # and the final_scf_energy is equal to the last value of scf_energies
        # the value of free_energy if list by ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):
        # which is the total DFT energy plus the Electronic entropic energy(very small).
        begin = self.run_info["scf_head_line"] + 2
        self.run_info["scf_energies"] = []
        for i in range(self.run_info["scf_steps"]):
            self.run_info["scf_energies"].append(float(self.lines[begin+i].split()[5]))
        # ----------------------------------------------------------------------
        # get the xyz structure from information extracted above:
        # WARNING: in low level print of cp2k, there is no structure coordinates
        # in the output file
        # ----------------------------------------------------------------------

    def plot_info(self):
        # now output the scf information in the current directory(post-processing)
        plt.plot(self.run_info["scf_energies"], marker="o")
        plt.title("Total energy per scf step")
        plt.xlabel("Scf step")
        plt.ylabel("Total energy(a.u.")
        plt.tight_layout()
        plt.savefig("total-energy-per-scf-step.png")
        plt.close()

    def markdown_report(self, md="scf-info.md"):
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# SCF实验统计\n")
            fout.write("## SCF参数\n")
            for item in self.scf_params:
                fout.write("- %s: %s\n" % (item, str(self.scf_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            # note that the accuracy for the seconds is not fully guranteed
            # e.g. 2019-11-26 12:09:36.487 is read as 2019-11-26 12:09:36
            start = datetime.datetime.strptime(self.run_info["start_time"].split()[7]+"-"+self.run_info["start_time"].split()[8].split(".")[0], "%Y-%m-%d-%H:%M:%S")
            stop = datetime.datetime.strptime(self.run_info["stop_time"].split()[7]+"-"+self.run_info["stop_time"].split()[8].split(".")[0], "%Y-%m-%d-%H:%M:%S")
            delta_t = stop -start
            fout.write("- Time consuming:\n")
            fout.write("  - totally %.1f seconds, or %.3f minutes or %.5f hours\n" % (delta_t.total_seconds(), delta_t.total_seconds()/60, delta_t.total_seconds()/3600))
            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            fout.write("![dft-energy-per-scf](total-energy-per-scf-step.png)")
    
    def to_json_string(self):
        """
        Note:
            self.run_info is a dict, but it cannot be dumped to json directory by json.dumps
            because there are objects inside self.run_info, like datetime that cannot
            be serialized by json.dumps()
        :return out: json string processed by json.dumps()
        """
        run_info = copy.deepcopy(self.run_info)
        scf_params = copy.deepcopy(self.scf_params)
        # convert datetime to str
        run_info["start_time"] = str(run_info["start_time"])
        # merge opt_params and run_info and output to json and return
        out = scf_params
        out.update(run_info)
        return json.dumps(out, indent=4)

    def export(self, directory):
        """
        """
        os.chdir(directory)
        os.system("mkdir -p post-processing")
        os.chdir("post-processing")
        self.plot_info()
        self.markdown_report()
        with open("scf.json", 'w') as fout:
            fout.write(self.to_json_string())
        os.chdir("../../")