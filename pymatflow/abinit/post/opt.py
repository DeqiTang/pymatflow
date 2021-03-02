#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import copy
import json
import datetime
import subprocess
import matplotlib.pyplot as plt

from pymatflow.base.atom import Atom


class Opt:
    """
    Note:
    """
    def __init__(self):
        """
        self.info: the dict containg all the information.
        self.lines: the opt output file as one string object, all information is parsed from it
        self.file: the abspath of the output  file of optimization run
        self.info["cells"]:
            a list of cell for every structure in self.trajectory.
            dimension: len(self.trajectory) * 9
        Note:
            opt will deal with the output of geometric optimization.
            it will automatically judge whether cell is optimized and build
            the cooresponding output.

            and even when the calculation is running or it is interrupted by
            you, you can post-process it!

            and only sigle dataset mode is supported
        """
        self.file = None
        self.info = {}
        self.info["cells"] = None #  optimized cell for every structure in trajectory
        self.info["acells"] = None #
        self.info["rprimds"] = None #
        self.trajectory = None


    def parse(self, file):
        """
        get the general information of opt run from opt run output file
        which is now stored in self.lines
        Note:
            self.get_outvars_before_and_after() will get the variables like
            typat and znucl and self.get_pseudo_info() will get the element
            name for every type of atom in typat, by analysing the pseudopotential
            file name for each type of atom.
            and self.get_trajectory() will use the information from self.outvars_before["typat"]
            and self.atom_type(by self.get_pseudo_info to set the atom name for every atom
            in every image in self.trajectory

            we don't read the final structure individually, because it is actually the final
            structure in self.trajectory, namely self.trajectory[-1]
        """
        self.file = os.path.abspath(file)

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()

        self.info = {} # make it a new empty dict, so that we can use get_info to get infor a brand new output
        self._get_outvars_before_and_after()
        self._get_pseduo_info() # also get information about the type of atom and its name
        self._get_run_info()
        self._get_trajectory()
        return self.info

    def _get_outvars_before_and_after(self):
        """
        """
        self.info["outvars"] = {}
        self.info["outvars"]["before"] = {}
        self.info["outvars"]["after"] = {}

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
                self.info["outvars"]["before"]["acell"] = [float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])]
            if self.lines[i].split()[0] == "diemac":
                self.info["outvars"]["before"]["diemac"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ecut":
                self.info["outvars"]["before"]["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ionmov":
                self.info["outvars"]["before"]["ionmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.info["outvars"]["before"]["istwfk"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "optcell":
                self.info["outvars"]["before"]["optcell"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "rprim":
                self.info["outvars"]["before"]["rprim"] = []
                for n in range(3):
                    self.info["outvars"]["before"]["rprim"].append(float(self.lines[i].split()[n+1]))
                for m in range(2):
                    for n in range(3):
                        self.info["outvars"]["before"]["rprim"].append(float(self.lines[i+m+1].split()[n]))
            if self.lines[i].split()[0] == "spgroup":
                self.info["outvars"]["before"]["spgroup"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "typat":
                self.info["outvars"]["before"]["typat"] = []
                for k in range(len(self.lines[i].split()) - 1):
                    self.info["outvars"]["before"]["typat"].append(int(self.lines[i].split()[k+1]))
                j = i + 1
                while self.lines[j].split()[0].isdigit(): # and len(self.lines[j].split()[1]) < 3:
                    for k in range(len(self.lines[j].split())):
                        self.info["outvars"]["before"]["typat"].append(int(self.lines[j].split()[k]))
                    j = j + 1
                #
            if self.lines[i].split()[0] == "znucl":
                self.info["outvars"]["before"]["znucl"] = []
                for k in range(len(self.lines[i].split()) - 1):
                    self.info["outvars"]["before"]["znucl"].append(float(self.lines[i].split()[k+1]))
                j = i + 1
                while len(self.lines[j].split()) != 0:
                    for k in range(len(self.lines[j].split())):
                        self.info["outvars"]["before"]["znucl"].append(float(self.lines[j].split()[k]))
                    j = j + 1
                #
        #
        for i in range(outvars_after_computation_line, len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0][0:3] == "===":
                break # end the final outvars after computation, end the loop.
            if self.lines[i].split()[0] == "acell":
                self.info["outvars"]["after"]["acell"] = [float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])]
            if self.lines[i].split()[0] == "diemac":
                self.info["outvars"]["after"]["diemac"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ecut":
                self.info["outvars"]["after"]["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ionmov":
                self.info["outvars"]["after"]["ionmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.info["outvars"]["after"]["istwfk"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "rprim":
                self.info["outvars"]["after"]["rprim"] = []
                for n in range(3):
                    self.info["outvars"]["after"]["rprim"].append(float(self.lines[i].split()[n+1]))
                for m in range(2):
                    for n in range(3):
                        self.info["outvars"]["after"]["rprim"].append(float(self.lines[i+m+1].split()[n]))
            if self.lines[i].split()[0] == "spgroup":
                self.info["outvars"]["after"]["spgroup"] = self.lines[i].split()[1]
            #
        #
        print("==============================================================\n")
        print("(by abinit.post.opt.opt_post.get_outvars_before_and_after()):\n")
        print("typat info:\n")
        print(self.info["outvars"]["before"]["typat"])
        print("you can check whether they are correct\n")

    def _get_pseduo_info(self):
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
        print("(by abinit.post.opt.opt_post.get_pseudo_info()):\n")
        print("Atom type info:\n")
        print(self.atom_type)
        print("you can check whether they are correct\n")

        #


    def _get_run_info(self):
        """
        self.info["iterations"]: scf iterations per ion step
        self.info["total_energies"]: total energies of every ion step
        self.info["fermi_energies"]: fermi energies of every ion step
        self.info["total_forces"]: total forces of every scf step
        """
        self.info["iterations"] = []
        self.info["total_energies"] = []

        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == "At" and self.lines[i].split()[1] == "SCF" and self.lines[i].split()[2] == "step":
                self.info["iterations"].append(int(self.lines[i].split()[3].split(",")[0]))
            if self.lines[i].split()[0] == "Total" and self.lines[i].split()[1] == "energy" and self.lines[i].split()[2] == "(etotal)":
                self.info["total_energies"].append(float(self.lines[i].split()[4])) # in unit of Hartree
            # get time information
            if self.lines[i].split()[0] == ".Starting" and self.lines[i].split()[1] == "date":
                self.info["start_time"] = datetime.datetime.strptime(self.lines[i].split()[4]+"-"+self.lines[i].split()[5]+"-"+self.lines[i].split()[6].split(".")[0]+"-"+self.lines[i+1].split()[3], "%d-%b-%Y-%Hh%M")
            # stop time is not available in output
            #if self.lines[i].split()[0] == ""
                #self.info["stop-time"] = datetime.datetime.strptime(self.lines[i].split()[4]+"-"+self.lines[i].split()[5]+"-"+self.lines[i].split()[6].split(".")[0]+"-"+self.lines[i+1].split()[3], "%d-%b-%Y-%Hh%M")
            #----------------------------------------------------------------------------------
            # note variable like ecut and ionmov and istwfk can appear twice in the output file
            # one in before the simulation, one after the simulation
            if self.lines[i].split()[0] == ".Version":
                self.info["version"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ecut":
                self.info["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ionmov":
                self.info["ionmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.info["istwfk"] = self.lines[i].split()[1]

    def _get_trajectory(self):
        """
        1 Bohr=0.5291772108 Angstroms
        xcart is in unit of bohr
        so we convert it to angstrom
        Note:
        """
        bohr = 0.5291772108
        self.info["trajectory"] = []
        self.info["cells"] = []
        self.info["acells"] = []
        self.info["rprimds"] = []
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            # used to get the trajectory
            if self.lines[i].split()[0] == "Cartesian" and self.lines[i].split()[1] == "coordinates" and self.lines[i].split()[2] == "(xcart)":
                atm = []
                j = i + 1
                while self.lines[j].split()[1] != "coordinates":
                    atm.append(Atom("XXX", float(self.lines[j].split()[0])*bohr, float(self.lines[j].split()[1])*bohr, float(self.lines[j].split()[2])*bohr))
                    j = j + 1
                self.info["trajectory"].append(atm)
            # --------------------------------
            # get the self.info["acells"] and self.info["rprimds"]
            # if not opt cell, the following two if code will never run to get the cell
            # so we can finally judget whether it is optcell by chekcing the length of self.info["rprimds"] and self.info["acells"]
            if self.lines[i].split()[0] == "Scale" and self.lines[i].split()[4] == "(acell)":
                self.info["acells"].append([float(self.lines[i+1].split()[0]), float(self.lines[i+1].split()[1]), float(self.lines[i+1].split()[2])])
            if self.lines[i].split()[0] == "Real" and self.lines[i].split()[1] == "space" and self.lines[i].split()[4] == "(rprimd)":
                rprimd = []
                for k in range(3):
                    for j in range(3):
                        rprimd.append(float(self.lines[i+1+k].split()[j]))
                self.info["rprimds"].append(rprimd)
            # end extract the cell
            # --------------------------------

        if len(self.info["acells"]) == len(self.info["trajectory"]):
            # so cell is optimized. acell and rprimd is print every ion step
            # Note rprimd is not rprim!!! and rprimd is actually rprim already scaled by acell
            # so to build cell from rprimd, we only need to multiply by bohr, and don't need to scale by acell
            self.info["cells"] = copy.deepcopy(self.info["rprimds"])
            for m in range(len(self.info["cells"])):
                for n in range(9):
                    self.info["cells"][m][n]= self.info["cells"][m][n] * bohr
            #

        # now using self.atom_type from self.get_pseudo_info()
        for i in range(len(self.info["trajectory"])):
            for j in range(len(self.info["trajectory"][i])):
                name = self.atom_type[self.info["outvars"]["before"]["typat"][j]]
                self.info["trajectory"][i][j].set_name(name)
        #

    def print_trajectory(self, directory="./"):
        with open(os.path.join(directory, "trajectory.xyz"), 'w') as fout:
            for i in range(len(self.info["trajectory"])):
                fout.write("%d\n" % len(self.info["trajectory"][i]))
                fout.write("i = %d\n" % i)
                for atom in self.info["trajectory"][i]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    def view_trajectory(self, trajfile="trajectory.xyz"):
        #os.system("xcrysden --xyz %s" % trajfile)
        subprocess.call(["xcrysden", "--xyz", trajfile])

    def print_final_structure(self, directory="./"):
        """
        we can finally judget whether it is optcell by chekcing the length of self.info["rprimds"] and self.info["acells"]
        against length of self.trajectory. if cell is optimized, in every iteration step of the optimization,
        it will print out the acell and rprimd and will be extract by self._get_trajectory.
        """
        bohr = 0.5291772108
        if len(self.info["acells"]) == 0 and len(self.info["rprimds"]) == 0:
            # not optimized cell
            cell = copy.deepcopy(self.info["outvars"]["before"]["rprim"])
            for n in range(3):
                cell[n] = cell[n] * self.info["outvars"]["before"]["acell"][0] * bohr
                cell[n+3] = cell[n+3] * self.info["outvars"]["before"]["acell"][1] * bohr
                cell[n+6] = cell[n+6] * self.info["outvars"]["before"]["acell"][2] * bohr
            #
        else:
            cell = self.info["cells"][-1]
        #
        with open(os.path.join(directory, "final-structure.xyz"), 'w') as fout:
            fout.write("%d\n" % len(self.info["trajectory"][-1]))
            fout.write("cell: %.9f %.9f %.9f | %.9f %.9f %.9f | %.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8]))
            for atom in self.info["trajectory"][-1]:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    def plot_info(self, directory="./"):
        """
        """
        plt.plot(self.info["iterations"])
        plt.title("Iterations per SCF")
        plt.xlabel("Scf cycles ")
        plt.ylabel("iterations")
        plt.tight_layout()
        plt.savefig(os.path.join(directory, "iterations-per-scf.png"))
        plt.close()

        plt.plot(self.info["total_energies"])
        plt.title("Total energies per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Total Energies (Hartree)")
        plt.tight_layout()
        plt.savefig(os.path.join(directory, "total_energies-per-scf.png"))
        plt.close()


    def markdown_report(self, directory="./"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(os.path.join(directory, "opt-info.md"), 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("## 几何优化参数\n")
            for item in self.info:
                fout.write("- %s: %s\n" % (item, str(self.info[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out

            # end the time information
            for item in self.info:
                fout.write("- %s: %s\n" % (item, str(self.info[item])))

            fout.write("## 运行信息图示\n")

            fout.write("Iterations per SCF\n")
            fout.write("![Iterations per SCF](iterations-per-scf.png)\n")

            fout.write("Total energies per scf\n")
            fout.write("![Total energies per SCF](total_energies-per-scf.png)\n")

    def to_json_string(self):
        """
        Note:
            self.info is a dict, but it cannot be dumped to json directory by json.dumps
            because there are objects inside self.info, like datetime and Atom() that cannot
            be serialized by json.dumps()
        :return out: json string processed by json.dumps()
        """
        out = copy.deepcopy(self.info)
        # convert Atoms to list object
        for image in out["trajectory"]:
            for i in range(len(image)):
                image[i] = [image[i].name, image[i].x, image[i].y, image[i].z]
        #
        # convert datetime to str
        out["start_time"] = str(out["start_time"])
        return json.dumps(out, indent=4)


    def export(self, directory="tmp-abinit-opt"):
        os.system("mkdir -p %s" % os.path.join(directory, "post-processing"))
        with open(os.path.join(directory, "post-processing/opt.json"), 'w') as fout:
            fout.write(self.to_json_string())
        self.print_trajectory(directory=os.path.join(directory, "post-processing"))
        self.print_final_structure(directory=os.path.join(directory, "post-processing"))
        self.plot_info(directory=os.path.join(directory, "post-processing"))
        self.markdown_report(directory=os.path.join(directory, "post-processing"))
