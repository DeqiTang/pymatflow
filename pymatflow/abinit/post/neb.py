#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import datetime
import subprocess
import matplotlib.pyplot as plt

from pymatflow.base.atom import Atom


class neb_out:
    """
    Note:
    """
    def __init__(self):
        """
        self.file:
            the output file of neb run.
        """
        self.file = None
        self.cell = None #
        self.neb_params = {}
        self.run_info = {}
        self.trajectory_initial = None
        self.trajectory_final = None



    def get_info(self, file):
        """
        get the general information of neb run from neb run output file
        which is now stored in self.lines
        """
        self.file = file
        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()

        self.get_neb_params_and_run_info()
        self.get_trajectory()
        return


    def get_trajectory(self):
        #
        outvars_after_computation_line = 0
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "after":
                outvars_after_computation_line = i
        #
        self.trajectory_final = []
        for i in range(outvars_after_computation_line, len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0].split("_")[0] == "xangst":
                atm = []
                # doesn't know name now
                atm.append(Atom("XXX", float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])))
                j = i + 1
                while len(self.lines[j].split()) == 3:
                    atm.append(Atom("XXX", float(self.lines[j].split()[0]), float(self.lines[j].split()[1]), float(self.lines[j].split()[2])))
                    j = j + 1
                self.trajectory_final.append(atm)

        #
        outvars_before_computation_line = 0
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "input":
                outvars_before_computation_line = i
        #
        self.trajectory_initial = []
        for i in range(outvars_before_computation_line, len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0].split("_")[0] == "xangst":
                atm = []
                # doesn't know name now
                atm.append(Atom("XXX", float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])))
                j = i + 1
                while len(self.lines[j].split()) == 3:
                    atm.append(Atom("XXX", float(self.lines[j].split()[0]), float(self.lines[j].split()[1]), float(self.lines[j].split()[2])))
                    j = j + 1
                self.trajectory_initial.append(atm)
        #

    def get_neb_params_and_run_info(self):
        """
        run_info["etotal-per-image"]: etotal of every image
        """
        self.run_info["etotal-per-image"] = []


        #
        outvars_after_computation_line = 0
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "after":
                outvars_after_computation_line = i
        #
        for i in range(outvars_after_computation_line, len(self.lines)):
            # if it is an empty line continue to next line
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == "ecut":
                self.neb_params["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0].split("_")[0] == "etotal":
                self.run_info["etotal-per-image"].append([self.lines[i].split()[0], float(self.lines[i].split()[1])])
            if self.lines[i].split()[0] == "imgmov":
                self.neb_params["imgmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.neb_params["istwfk"] = self.lines[i].split()[1]

        # get time information
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == ".Starting" and self.lines[i].split()[1] == "date":
                self.run_info["start_time"] = datetime.datetime.strptime(self.lines[i].split()[4]+"-"+self.lines[i].split()[5]+"-"+self.lines[i].split()[6].split(".")[0]+"-"+self.lines[i+1].split()[3], "%d-%b-%Y-%Hh%M")
            # stop time is not available in output
            #if self.lines[i].split()[0] == ""
                #self.run_info["stop-time"] = datetime.datetime.strptime(self.lines[i].split()[4]+"-"+self.lines[i].split()[5]+"-"+self.lines[i].split()[6].split(".")[0]+"-"+self.lines[i+1].split()[3], "%d-%b-%Y-%Hh%M")



class neb_post:
    """
    Note:
    """
    def __init__(self, output):
        """
        output:
            the output file of neb run.
        """
        self.file = output
        self.cell = None #
        self.neb_params = {}
        self.run_info = {}
        self.trajectory_initial = None
        self.trajectory_final = None

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_info()

    def get_info(self):
        """
        get the general information of neb run from neb run output file
        which is now stored in self.lines
        """

        self.get_neb_params_and_run_info()
        self.get_trajectory()
        return


    def get_trajectory(self):
        #
        outvars_after_computation_line = 0
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "after":
                outvars_after_computation_line = i
        #
        self.trajectory_final = []
        for i in range(outvars_after_computation_line, len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0].split("_")[0] == "xangst":
                atm = []
                # doesn't know name now
                atm.append(Atom("XXX", float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])))
                j = i + 1
                while len(self.lines[j].split()) == 3:
                    atm.append(Atom("XXX", float(self.lines[j].split()[0]), float(self.lines[j].split()[1]), float(self.lines[j].split()[2])))
                    j = j + 1
                self.trajectory_final.append(atm)

        #
        outvars_before_computation_line = 0
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "input":
                outvars_before_computation_line = i
        #
        self.trajectory_initial = []
        for i in range(outvars_before_computation_line, len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0].split("_")[0] == "xangst":
                atm = []
                # doesn't know name now
                atm.append(Atom("XXX", float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])))
                j = i + 1
                while len(self.lines[j].split()) == 3:
                    atm.append(Atom("XXX", float(self.lines[j].split()[0]), float(self.lines[j].split()[1]), float(self.lines[j].split()[2])))
                    j = j + 1
                self.trajectory_initial.append(atm)
        #

    def get_neb_params_and_run_info(self):
        """
        run_info["etotal-per-image"]: etotal of every image
        """
        self.run_info["etotal-per-image"] = []


        #
        outvars_after_computation_line = 0
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "after":
                outvars_after_computation_line = i
        #
        for i in range(outvars_after_computation_line, len(self.lines)):
            # if it is an empty line continue to next line
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == "ecut":
                self.neb_params["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0].split("_")[0] == "etotal":
                self.run_info["etotal-per-image"].append([self.lines[i].split()[0], float(self.lines[i].split()[1])])
            if self.lines[i].split()[0] == "imgmov":
                self.neb_params["imgmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.neb_params["istwfk"] = self.lines[i].split()[1]

        # get time information
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == ".Starting" and self.lines[i].split()[1] == "date":
                self.run_info["start-time"] = datetime.datetime.strptime(self.lines[i].split()[4]+"-"+self.lines[i].split()[5]+"-"+self.lines[i].split()[6].split(".")[0]+"-"+self.lines[i+1].split()[3], "%d-%b-%Y-%Hh%M")
            # stop time is not available in output
            #if self.lines[i].split()[0] == ""
                #self.run_info["stop-time"] = datetime.datetime.strptime(self.lines[i].split()[4]+"-"+self.lines[i].split()[5]+"-"+self.lines[i].split()[6].split(".")[0]+"-"+self.lines[i+1].split()[3], "%d-%b-%Y-%Hh%M")


    def print_trajectory(self, xyz_initial="trajectory-initial.xyz", xyz_final="trajectory-final.xyz"):
        with open(xyz_initial, 'w') as fout:
            for i in range(len(self.trajectory_initial)):
                fout.write("%d\n" % len(self.trajectory_initial[i]))
                fout.write("i = %d\n" % i)
                for atom in self.trajectory_initial[i]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
        with open(xyz_final, 'w') as fout:
            for i in range(len(self.trajectory_final)):
                fout.write("%d\n" % len(self.trajectory_final[i]))
                fout.write("i = %d\n" % i)
                for atom in self.trajectory_final[i]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    def view_trajectory(self, trajfile_initial="trajectory-initial.xyz", trajfile_final="trajectory-final.xyz"):
        #os.system("xcrysden --xyz %s" % trajfile)
        subprocess.call(["xcrysden", "--xyz", trajfile_initial])
        subprocess.call(["xcrysden", "--xyz", trajfile_final])

    def plot_run_info(self):
        """
        """

        plt.plot([self.run_info["etotal-per-image"][i][1] for i in range(len(self.run_info["etotal-per-image"]))])
        plt.title("Total energies per image")
        plt.xlabel("image")
        plt.ylabel("Total Energies (Hartree)")
        plt.tight_layout()
        plt.savefig("etotal-per-image.png")
        plt.close()


    def markdown_report(self, md="TransitionStateSearchReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 过渡态搜索优化实验统计\n")
            fout.write("## 过渡态参数\n")
            for item in self.neb_params:
                fout.write("- %s: %s\n" % (item, str(self.neb_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out

            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")

            fout.write("Total energies per image\n")
            fout.write("![Total energies per image](etotal-per-image.png)\n")


    def export(self):
        #self.to_xyz()
        self.print_trajectory()
        self.plot_run_info()
        self.markdown_report("TransitionStateSearchReport.md")
