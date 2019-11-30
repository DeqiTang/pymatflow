#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import datetime
import subprocess
import matplotlib.pyplot as plt

from pymatflow.base.atom import Atom


class opt_post:
    """
    Note:
        opt_post can extract information for the geometric optimization running,
        including 'relax' and 'vc-relax'. 
        it will printout the trajectory file(xyz format) and the final optimized
        structure(if not relaxed, the final structure of the running).
        so even when your ion step not converged within maximum steps, it can also
        extract the structure in the final ion step and print it out.
    """
    def __init__(self, output, run_type):
        """
        output:
            the output file of opt run. it could be a geo opt or a cell opt
        
        relaxed:
            whether the structure successfully relaxed.
        """
        self.file = output
        self.run_type = run_type 
        self.cell = None # optimized cell
        self.atoms = None  # optimized atoms
        self.opt_params = {}
        self.run_info = {}
        self.trajectory = None
        self.relaxed = None # whether structure is relaxed or vcr-relaxed

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_info()
        self.get_trajectory()

    def get_info(self):
        """
        get the general information of opt run from opt run output file
        which is now stored in self.lines
        """
        # check whether successfully relaxed
        self.relaxed = False
        for line in self.lines:
            if line == "Begin final coordinates\n":
                self.relaxed = True
                break
        #
        if self.relaxed == True:
            if self.run_type == "relax":
                self.get_structure_relax()
            elif self.run_type == "vc-relax":
                self.get_structure_vc_relax()
        
        self.get_opt_params_and_run_info()

    def get_structure_relax(self):
        """
        """
        self.atoms = []
        # get the line number of the 'Begin final coordinates'
        # and 'End final coordinates'
        begin_final_coord_line = 0
        end_final_coord_line = 0
        while self.lines[begin_final_coord_line] != "Begin final coordinates\n":
            begin_final_coord_line += 1
        while self.lines[end_final_coord_line] != "End final coordinates\n":
            end_final_coord_line += 1

        # coords(in relax running it will not print the cell as it does not change)
        for i in range(begin_final_coord_line+3, end_final_coord_line):
            self.atoms.append(Atom(self.lines[i].split()[0], float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])))
 
    def get_structure_vc_relax(self):
        """
        """
        self.cell = []
        self.atoms = []
        # get the line number of the 'Begin final coordinates'
        # and 'End final coordinates'
        begin_final_coord_line = 0
        end_final_coord_line = 0
        while self.lines[begin_final_coord_line] != "Begin final coordinates\n":
            begin_final_coord_line += 1
        while self.lines[end_final_coord_line] != "End final coordinates\n":
            end_final_coord_line += 1

        # get cell and coords
        self.cell = []
        for i in range(begin_final_coord_line+5, begin_final_coord_line+8):
            for j in range(3):
                self.cell.append(float(self.lines[i].split()[j]))
        for i in range(begin_final_coord_line+10, end_final_coord_line):
            self.atoms.append(Atom(self.lines[i].split()[0], float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])))
   

    #

    def get_trajectory(self):
        self.trajectory = []
        for i in range(len(self.lines)):
            if len(self.lines[i].split()) > 0 and self.lines[i].split()[0] == "ATOMIC_POSITIONS":
                atm = []
                j = i + 1
                while len(self.lines[j].split()) == 4:
                    atm.append(Atom(self.lines[j].split()[0], float(self.lines[j].split()[1]), float(self.lines[j].split()[2]), float(self.lines[j].split()[3])))
                    j = j + 1
                self.trajectory.append(atm)
                

    def get_opt_params_and_run_info(self):
        """
        run_info["iterations"]: scf iterations per scf step
        run_info["total-energies"]: total energies of every scf step
        run_info["fermi-energies"]: fermi energies of every scf step
        run_info["total-forces"]: total forces of every scf step
        """
        self.run_info["iterations"] = []
        self.run_info["total-energies"] = []
        self.run_info["fermi-energies"] = []
        self.run_info["total-forces"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
            if line.split()[0] == "Program" and line.split()[1] == "PWSCF" and line.split()[3] == "starts":
                self.run_info["start-time"] = line.split("\n")[0]
            if line.split()[0] == "This" and line.split()[1] == "run" and line.split()[3] == "terminated":
                self.run_info["stop-time"] = line.split("\n")[0]
            if line.split()[0] == "kinetic-energy":
                self.opt_params["ecutwfc"] = int(float(line.split()[3]))
            if line.split()[0] == "convergence threshold":
                self.opt_params["conv_thr"] = float(line.split()[3])
            if line.split()[0] == "mixing" and line.split()[1] == 'beta':
                self.opt_params["mixing_beta"] = float(line.split()[3])
            if line.split()[0] == "number" and line.split()[2] == 'k':
                self.opt_params["degauss"] = float(line.split()[9])
            if line.split()[0] == "convergence" and line.split()[3] == "achieved":
                self.run_info["iterations"].append(int(line.split()[5]))
            if line.split()[0] == "!" and line.split()[5] == "Ry":
                self.run_info["total-energies"].append(float(line.split()[4]))
            if line.split()[0] ==  "the" and line.split()[1] == "Fermi":
                self.run_info["fermi-energies"].append(float(line.split()[4]))
            if line.split()[0] == "Total" and line.split()[1] == "force":
                self.run_info["total-forces"].append(float(line.split()[3]))

        self.run_info["scf-cycles"] = len(self.run_info["iterations"])
        if self.run_type == "relax":
            self.run_info["ion-steps"] = len(self.run_info["iterations"]) - 1
        elif self.run_type == "vc-relax":
            self.run_info["ion-steps"] = len(self.run_info["iterations"]) - 2

    def to_xyz(self, xyz="optimized.xyz"):
        if self.relaxed == False:
            with open("final-structure(not-relaxed).xyz", 'w') as fout:
                fout.write("%d\n" % len(self.trajectory[0]))
                fout.write("Warning(%s): structure failed to be relaxed or vc-relaxed, this is the final structure(unrelaxed)\n" % self.run_type)
                for atom in self.trajectory[-1]:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
            return
        # printout relaxed structure
        cell = self.cell
        with open(xyz, 'w') as fout:
            fout.write("%d\n" % len(self.atoms))
            if self.run_type == "vc-relax":
                fout.write("cell: %.9f %.9f %.9f | %.9f %.9f %.9f | %.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8]))
            else:
                fout.write("type of opt run: relax -> the cell is not changed, so go and find the original cell\n")
            for atom in self.atoms:
                fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))

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

    def plot_run_info(self):
        """
        """
        plt.plot(self.run_info["iterations"])
        plt.title("Iterations per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("iterations")
        plt.tight_layout()
        plt.savefig("iterations-per-scf.png")
        plt.close()

        plt.plot(self.run_info["total-energies"])
        plt.title("Total energies per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Total Energies (Ry)")
        plt.tight_layout()
        plt.savefig("total-energies-per-scf.png")
        plt.close()

        plt.plot(self.run_info["fermi-energies"])
        plt.title("Fermi energies per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Fermi energies (eV)")
        plt.tight_layout()
        plt.savefig("fermi-energies-per-scf.png")
        plt.close()

        plt.plot(self.run_info["total-forces"])
        plt.title("Total forces per SCF")
        plt.xlabel("Scf cycles")
        plt.ylabel("Total forces (Ry/au)")
        plt.tight_layout()
        plt.savefig("total-forces-per-scf.png")
        plt.close()

    def markdown_report(self, md="OptimizationReport.md"):
        """
        when writing Chinese to a file you must specify
        encoding='utf-8' when open the file for writing
        """
        with open(md, 'w', encoding='utf-8') as fout:
            fout.write("# 几何优化实验统计\n")
            fout.write("几何优化类型: %s\n" % self.run_type)
            fout.write("是否成功优化: %s\n" % str(self.relaxed))
            fout.write("## 优化参数\n")
            for item in self.opt_params:
                fout.write("- %s: %s\n" % (item, str(self.opt_params[item])))
            fout.write("## 运行信息\n")
            # calculate the running time and print it out
            # Importante: the length of the time string might be different, depending
            # on the value of hours and minutes and seconds. if they are two digits
            # number, they will be divided like: '11: 6: 2', only when they all are
            # two digtis number, they will not be divided '11:16:12'
            # so we have to preprocess it to build the right time string to pass into 
            # datetime.datetime.strptime()
            if len(self.run_info["start-time"].split()) == 8:
                start_str = self.run_info["start-time"].split()[5]+"-"+self.run_info["start-time"].split()[7]
            elif len(self.run_info["start-time"].split()) == 9:
                start_str = self.run_info["start-time"].split()[5]+"-"+self.run_info["start-time"].split()[7]+self.run_info["start-time"].split()[8]
            elif len(self.run_info["start-time"].split()) == 10:
                start_str = self.run_info["start-time"].split()[5]+"-"+self.run_info["start-time"].split()[7]+self.run_info["start-time"].split()[8]+self.run_info["start-time"].split()[9]
            else:
                print("===============================================\n")
                print("                  Warning !!!\n")
                print("===============================================\n")
                print("qe.post.opt.markdown_report:\n")
                print("failed to parse start-time string\n")
                sys.exit(1)
            if len(self.run_info["stop-time"].split()) == 7:
                stop_str = self.run_info["stop-time"].split()[6]+"-"+self.run_info["stop-time"].split()[5]
            elif len(self.run_info["stop-time"].split()) == 8:
                stop_str = self.run_info["stop-time"].split()[7]+"-"+self.run_info["stop-time"].split()[5]+self.run_info["stop-time"].split()[6]
            elif len(self.run_info["stop-time"].split()) == 9:
                stop_str = self.run_info["stop-time"].split()[8]+"-"+self.run_info["stop-time"].split()[5]+self.run_info["stop-time"].split()[6]+self.run_info["stop-time"].split()[7]
            else:
                print("===============================================\n")
                print("                  Warning !!!\n")
                print("===============================================\n")
                print("qe.post.opt.markdown_report:\n")
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
            fout.write("Iterations per SCF\n")
            fout.write("![Iterations per SCF](iterations-per-scf.png)\n")
            
            fout.write("Total energies per SCF\n")
            fout.write("![Total energies per SCF](total-energies-per-scf.png)\n")

            fout.write("Fermi energies per SCF\n")
            fout.write("![Fermi energies per SCF](fermi-energies-per-scf.png)\n")

            fout.write("Total forces per SCF\n")
            fout.write("![Total forces per SCF](total-forces-per-scf.png)\n")


    def export(self):
        self.to_xyz()
        self.print_trajectory()
        self.plot_run_info()
        self.markdown_report("OptimizationReport.md")
