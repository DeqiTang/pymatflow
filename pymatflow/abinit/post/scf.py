
import os
import copy
import datetime
import subprocess
import matplotlib.pyplot as plt

from pymatflow.base.atom import Atom



class scf_out:
    """
    Note:
    """
    def __init__(self):
        """
        self.file:
            the output file of static run.
        """
        self.file = None
        self.scf_params = {}
        self.run_info = {}
        self.outvars_before = {}
        self.outvars_after = {}

    def get_info(self, file):
        """
        get the general information of opt run from static run output file
        which is now stored in self.lines
        """
        self.file = file
        with open(self.file, 'r') as fin:
            self.lines = fout.readlines()

        self.get_outvars_before_and_after()
        self.get_pseudo_info() # also get information about the type of atom and its name
        self.get_scf_params_and_run_info()
        return

    def get_outvars_before_and_after(self):
        """
        """
        self.outvars_before = {}
        self.outvars_after = {}
        outvars_before_computation_line = 0
        outvars_after_computation_line = 0
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "input":
                outvars_before_computation_line = i
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "after":
                outvars_after_computation_line = i
        #
        for i in range(outvars_before_computation_line, len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0][0:3] == "===":
                break # end the first outvars before computation, end the loop.
            if self.lines[i].split()[0] == "acell":
                self.outvars_before["acell"] = [float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])]
            if self.lines[i].split()[0] == "diemac":
                self.outvars_before["diemac"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ecut":
                self.outvars_before["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ionmov":
                self.outvars_before["ionmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.outvars_before["istwfk"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "optcell":
                self.outvars_before["optcell"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "rprim":
                self.outvars_before["rprim"] = []
                for n in range(3):
                    self.outvars_before["rprim"].append(float(self.lines[i].split()[n+1]))
                for m in range(2):
                    for n in range(3):
                        self.outvars_before["rprim"].append(float(self.lines[i+m+1].split()[n]))
            if self.lines[i].split()[0] == "spgroup":
                self.outvars_before["spgroup"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "typat":
                self.outvars_before["typat"] = []
                for k in range(len(self.lines[i].split()) - 1):
                    self.outvars_before["typat"].append(int(self.lines[i].split()[k+1]))
                j = i + 1
                while self.lines[j].split()[0] != "-":
                    for k in range(len(self.lines[j].split())):
                        self.outvars_before["typat"].append(int(self.lines[j].split()[k]))
                    j = j + 1
                #
            if self.lines[i].split()[0] == "znucl":
                self.outvars_before["znucl"] = []
                for k in range(len(self.lines[i].split()) - 1):
                    self.outvars_before["znucl"].append(float(self.lines[i].split()[k+1]))
                j = i + 1
                while len(self.lines[j].split()) != 0:
                    for k in range(len(self.lines[j].split())):
                        self.outvars_before["znucl"].append(float(self.lines[j].split()[k]))
                    j = j + 1
                #
        #
        for i in range(outvars_after_computation_line, len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0][0:3] == "===":
                break # end the final outvars after computation, end the loop.
            if self.lines[i].split()[0] == "acell":
                self.outvars_after["acell"] = [float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])]
            if self.lines[i].split()[0] == "diemac":
                self.outvars_after["diemac"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ecut":
                self.outvars_after["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ionmov":
                self.outvars_after["ionmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.outvars_after["istwfk"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "rprim":
                self.outvars_after["rprim"] = []
                for n in range(3):
                    self.outvars_after["rprim"].append(float(self.lines[i].split()[n+1]))
                for m in range(2):
                    for n in range(3):
                        self.outvars_after["rprim"].append(float(self.lines[i+m+1].split()[n]))
            if self.lines[i].split()[0] == "spgroup":
                self.outvars_after["spgroup"] = self.lines[i].split()[1]
            #
        #
        print("==============================================================\n")
        print("(by abinit.post.scf.scf_post.get_outvars_before_and_after()):\n")
        print("typat info:\n")
        print(self.outvars_before["typat"])
        print("you can check whether they are correct\n")

    def get_pseudo_info(self):
        """
        also get information about the type of atom and its name
        self.atom_type:
            {1: 'element1', 2: 'element2', ...}
        #
        """
        self.atom_type = {}
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0 or len(self.lines[i].split()) == 1:
                continue
            if self.lines[i].split()[0] == "-" and self.lines[i].split()[1] == "pspini:":
                self.atom_type[int(self.lines[i].split()[4])] = self.lines[i].split()[8].split(".")[0]

        print("=====================================================\n")
        print("(by abinit.post.scf.scf_post.get_pseudo_info()):\n")
        print("Atom type info:\n")
        print(self.atom_type)
        print("you can check whether they are correct\n")

        #


    def get_scf_params_and_run_info(self):
        """
        """
        #self.run_info["iterations"] = []
        #self.run_info["total-energies"] = []

        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            # get time information
            if self.lines[i].split()[0] == ".Starting" and self.lines[i].split()[1] == "date":
                self.run_info["start_time"] = datetime.datetime.strptime(self.lines[i].split()[4]+"-"+self.lines[i].split()[5]+"-"+self.lines[i].split()[6].split(".")[0]+"-"+self.lines[i+1].split()[3], "%d-%b-%Y-%Hh%M")
            # stop time is not available in output
            #if self.lines[i].split()[0] == ""
                #self.run_info["stop-time"] = datetime.datetime.strptime(self.lines[i].split()[4]+"-"+self.lines[i].split()[5]+"-"+self.lines[i].split()[6].split(".")[0]+"-"+self.lines[i+1].split()[3], "%d-%b-%Y-%Hh%M")
            #----------------------------------------------------------------------------------
            # note variable like ecut and ionmov and istwfk can appear twice in the output file
            # one in before the simulation, one after the simulation
            if self.lines[i].split()[0] == "ecut":
                self.scf_params["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ionmov":
                self.scf_params["ionmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.scf_params["istwfk"] = self.lines[i].split()[1]



class scf_post:
    """
    Note:
    """
    def __init__(self, output):
        """
        output:
            the output file of static run.
        """
        self.file = output
        self.scf_params = {}
        self.run_info = {}
        self.outvars_before = {}
        self.outvars_after = {}

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_info()

    def get_info(self):
        """
        get the general information of opt run from static run output file
        which is now stored in self.lines

        """
        self.get_outvars_before_and_after()
        self.get_pseduo_info() # also get information about the type of atom and its name
        self.get_scf_params_and_run_info()
        return

    def get_outvars_before_and_after(self):
        """
        """
        self.outvars_before = {}
        self.outvars_after = {}
        outvars_before_computation_line = 0
        outvars_after_computation_line = 0
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "input":
                outvars_before_computation_line = i
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "after":
                outvars_after_computation_line = i
        #
        for i in range(outvars_before_computation_line, len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0][0:3] == "===":
                break # end the first outvars before computation, end the loop.
            if self.lines[i].split()[0] == "acell":
                self.outvars_before["acell"] = [float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])]
            if self.lines[i].split()[0] == "diemac":
                self.outvars_before["diemac"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ecut":
                self.outvars_before["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ionmov":
                self.outvars_before["ionmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.outvars_before["istwfk"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "optcell":
                self.outvars_before["optcell"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "rprim":
                self.outvars_before["rprim"] = []
                for n in range(3):
                    self.outvars_before["rprim"].append(float(self.lines[i].split()[n+1]))
                for m in range(2):
                    for n in range(3):
                        self.outvars_before["rprim"].append(float(self.lines[i+m+1].split()[n]))
            if self.lines[i].split()[0] == "spgroup":
                self.outvars_before["spgroup"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "typat":
                self.outvars_before["typat"] = []
                for k in range(len(self.lines[i].split()) - 1):
                    self.outvars_before["typat"].append(int(self.lines[i].split()[k+1]))
                j = i + 1
                while self.lines[j].split()[0] != "-":
                    for k in range(len(self.lines[j].split())):
                        self.outvars_before["typat"].append(int(self.lines[j].split()[k]))
                    j = j + 1
                #
            if self.lines[i].split()[0] == "znucl":
                self.outvars_before["znucl"] = []
                for k in range(len(self.lines[i].split()) - 1):
                    self.outvars_before["znucl"].append(float(self.lines[i].split()[k+1]))
                j = i + 1
                while len(self.lines[j].split()) != 0:
                    for k in range(len(self.lines[j].split())):
                        self.outvars_before["znucl"].append(float(self.lines[j].split()[k]))
                    j = j + 1
                #
        #
        for i in range(outvars_after_computation_line, len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0][0:3] == "===":
                break # end the final outvars after computation, end the loop.
            if self.lines[i].split()[0] == "acell":
                self.outvars_after["acell"] = [float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])]
            if self.lines[i].split()[0] == "diemac":
                self.outvars_after["diemac"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ecut":
                self.outvars_after["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ionmov":
                self.outvars_after["ionmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.outvars_after["istwfk"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "rprim":
                self.outvars_after["rprim"] = []
                for n in range(3):
                    self.outvars_after["rprim"].append(float(self.lines[i].split()[n+1]))
                for m in range(2):
                    for n in range(3):
                        self.outvars_after["rprim"].append(float(self.lines[i+m+1].split()[n]))
            if self.lines[i].split()[0] == "spgroup":
                self.outvars_after["spgroup"] = self.lines[i].split()[1]
            #
        #
        print("==============================================================\n")
        print("(by abinit.post.scf.scf_post.get_outvars_before_and_after()):\n")
        print("typat info:\n")
        print(self.outvars_before["typat"])
        print("you can check whether they are correct\n")

    def get_pseduo_info(self):
        """
        also get information about the type of atom and its name
        self.atom_type:
            {1: 'element1', 2: 'element2', ...}
        #
        """
        self.atom_type = {}
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0 or len(self.lines[i].split()) == 1:
                continue
            if self.lines[i].split()[0] == "-" and self.lines[i].split()[1] == "pspini:":
                self.atom_type[int(self.lines[i].split()[4])] = self.lines[i].split()[8].split(".")[0]

        print("=====================================================\n")
        print("(by abinit.post.scf.scf_post.get_pseudo_info()):\n")
        print("Atom type info:\n")
        print(self.atom_type)
        print("you can check whether they are correct\n")

        #


    def get_scf_params_and_run_info(self):
        """
        """
        #self.run_info["iterations"] = []
        #self.run_info["total-energies"] = []

        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            # get time information
            if self.lines[i].split()[0] == ".Starting" and self.lines[i].split()[1] == "date":
                self.run_info["start-time"] = datetime.datetime.strptime(self.lines[i].split()[4]+"-"+self.lines[i].split()[5]+"-"+self.lines[i].split()[6].split(".")[0]+"-"+self.lines[i+1].split()[3], "%d-%b-%Y-%Hh%M")
            # stop time is not available in output
            #if self.lines[i].split()[0] == ""
                #self.run_info["stop-time"] = datetime.datetime.strptime(self.lines[i].split()[4]+"-"+self.lines[i].split()[5]+"-"+self.lines[i].split()[6].split(".")[0]+"-"+self.lines[i+1].split()[3], "%d-%b-%Y-%Hh%M")
            #----------------------------------------------------------------------------------
            # note variable like ecut and ionmov and istwfk can appear twice in the output file
            # one in before the simulation, one after the simulation
            if self.lines[i].split()[0] == "ecut":
                self.scf_params["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ionmov":
                self.scf_params["ionmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.scf_params["istwfk"] = self.lines[i].split()[1]


    def plot_run_info(self):
        """
        """
        #plt.plot(self.run_info["iterations"])
        #plt.title("Iterations per SCF")
        #plt.xlabel("Scf cycles ")
        #plt.ylabel("iterations")
        #plt.tight_layout()
        #plt.savefig("iterations-per-scf.png")
        #plt.close()
        pass


    def markdown_report(self, md="StaticCalcReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 静态计算实验统计\n")
            fout.write("## 计算参数\n")
            for item in self.scf_params:
                fout.write("- %s: %s\n" % (item, str(self.scf_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out

            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")

            fout.write("Iterations per SCF\n")
            fout.write("![Iterations per SCF](iterations-per-scf.png)\n")


    def export(self):
        self.plot_run_info()
        self.markdown_report("StaticCalcReport.md")
