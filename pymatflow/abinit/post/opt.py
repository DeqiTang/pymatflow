#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import datetime
import subprocess
import matplotlib.pyplot as plt

from pymatflow.base.atom import Atom


class opt_post:
    """
    Note:
    """
    def __init__(self, output):
        """
        output:
            the output file of optimization run.
        """
        self.file = output
        self.cell = None #  optimized cell
        self.structure_final = None
        self.opt_params = {}
        self.run_info = {}
        self.trajectory = None

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_info()

    def get_info(self):
        """
        get the general information of opt run from opt run output file
        which is now stored in self.lines
        """
        
        self.get_opt_params_and_run_info()
        self.get_trajectory()
        self.get_final_structure()
        return



    def get_opt_params_and_run_info(self):
        """
        run_info["iterations"]: scf iterations per scf step
        run_info["total-energies"]: total energies of every scf step
        run_info["fermi-energies"]: fermi energies of every scf step
        run_info["total-forces"]: total forces of every scf step
        """
        self.run_info["iterations"] = []
        self.run_info["total-energies"] = []

        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == "At" and self.lines[i].split()[1] == "SCF" and self.lines[i].split()[6] == "converged":
                self.run_info["iterations"].append(int(self.lines[i].split()[3].split(",")[0]))
            if self.lines[i].split()[0] == "Total" and self.lines[i].split()[1] == "energy" and self.lines[i].split()[2] == "(etotal)":
                self.run_info["total-energies"].append(float(self.lines[i].split()[4])) # in unit of Hartree
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
                self.opt_params["ecut"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "ionmov":
                self.opt_params["ionmov"] = self.lines[i].split()[1]
            if self.lines[i].split()[0] == "istwfk":
                self.opt_params["istwfk"] = self.lines[i].split()[1]

    def get_trajectory(self):
        """
        1 Bohr=0.5291772108 Angstroms
        xcart is in unit of bohr
        so we convert it to angstrom
        """
        bohr = 0.5291772108
        self.trajectory = []
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "Cartesian" and self.lines[i].split()[1] == "coordinates" and self.lines[i].split()[2] == "(xcart)":
                atm = []
                j = i + 1
                while self.lines[j].split()[1] != "coordinates":
                    atm.append(Atom("XXX", float(self.lines[j].split()[0])*bohr, float(self.lines[j].split()[1])*bohr, float(self.lines[j].split()[2])*bohr))
                    j = j + 1
                self.trajectory.append(atm)


    def get_final_structure(self):
        #
        outvars_after_computation_line = 0
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "-outvars:" and self.lines[i].split()[5] == "after":
                outvars_after_computation_line = i
        #
        for i in range(outvars_after_computation_line, len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            # get the final structure from outvars after computation
            if self.lines[i].split()[0] == "xangst":
                self.structure_final = []
                # doesn't know name now
                self.structure_final.append(Atom("XXX", float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])))
                j = i + 1
                while len(self.lines[j].split()) == 3:
                    self.structure_final.append(Atom("XXX", float(self.lines[j].split()[0]), float(self.lines[j].split()[1]), float(self.lines[j].split()[2])))
                    j = j + 1
                #
            # get the cell by acell and rprim
            if self.lines[i].split()[0] == "acell":
                self.run_info["acell-final"] = [float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])]
            if self.lines[i].split()[0] == "rprim":
                self.cell = []
                for j in range(3):
                    self.cell.append(float(self.lines[i].split()[j+1]))
                for j in range(3):
                    self.cell.append(float(self.lines[i+1].split()[j]))     
                for j in range(3):
                    self.cell.append(float(self.lines[i+2].split()[j]))
        # now we must multiply self.cell with run_info["acell-final"] and bohr to build the real angstrom cell
        bohr = 0.5291772108
        for k in range(3):
            self.cell[k] = self.cell[k] * self.run_info["acell-final"][0] * bohr
            self.cell[k+3] = self.cell[k+3] * self.run_info["acell-final"][1] * bohr
            self.cell[k+6] = self.cell[k+6] * self.run_info["acell-final"][2] * bohr
        # end extract the cell

        
   
    def print_trajectory(self, xyz="trajectory.xyz"):
        with open(xyz, 'w') as fout:
            for i in range(len(self.trajectory)):
                fout.write("%d\n" % len(self.trajectory[i]))
                fout.write("i = %d\n" % i)
                for atom in self.trajectory[i]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    def view_trajectory(self, trajfile="trajectory.xyz"):
        #os.system("xcrysden --xyz %s" % trajfile)
        subprocess.call(["xcrysden", "--xyz", trajfile])

    def print_final_structure(self, xyz="optimized.xyz"):
        #if self.relaxed == False:
        #    with open("final-structure(not-relaxed).xyz", 'w') as fout:
        #        fout.write("%d\n" % len(self.trajectory[0]))
        #        fout.write("Warning(%s): structure failed to be relaxed or vc-relaxed, this is the final structure(unrelaxed)\n" % self.run_type)
        #        for atom in self.trajectory[-1]:
        #            fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
        #    return

        # printout relaxed structure
        cell = self.cell
        with open(xyz, 'w') as fout:
            fout.write("%d\n" % len(self.structure_final))
            fout.write("cell: %.9f %.9f %.9f | %.9f %.9f %.9f | %.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8]))
            for atom in self.structure_final:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

    def plot_run_info(self):
        """
        """
        plt.plot(self.run_info["iterations"])
        plt.title("Iterations per SCF")
        plt.xlabel("Scf cycles ")
        plt.ylabel("iterations")
        plt.tight_layout()
        plt.savefig("iterations-per-scf.png")
        plt.close()

        plt.plot(self.run_info["total-energies"])
        plt.title("Total energies per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Total Energies (Hartree)")
        plt.tight_layout()
        plt.savefig("total-energies-per-scf.png")
        plt.close()


    def markdown_report(self, md="GeometricOptimizationReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("## 几何优化参数\n")
            for item in self.opt_params:
                fout.write("- %s: %s\n" % (item, str(self.opt_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out

            # end the time information
            for item in self.run_info:
                fout.write("- %s: %s\n" % (item, str(self.run_info[item])))

            fout.write("## 运行信息图示\n")
            
            fout.write("Iterations per SCF\n")
            fout.write("![Iterations per SCF](iterations-per-scf.png)\n")

            fout.write("Total energies per scf\n")
            fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")


    def export(self):
        self.print_trajectory()
        self.print_final_structure()
        self.plot_run_info()
        self.markdown_report("GeometricOptimizationReport.md")
