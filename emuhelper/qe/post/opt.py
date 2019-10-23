#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import matplotlib.pyplot as plt

from emuhelper.base.atom import Atom


class opt_post:
    """
    """
    def __init__(self, output, run_type):
        """
        output is the output file of opt run. it could be a geo opt
        or a cell opt
        """
        self.file = output
        self.run_type = run_type 
        self.cell = None # optimized cell
        self.atoms = None  # optimized atoms
        self.opt_params = {}
        self.run_info = {}

        with open(self.file, 'r') as fout:
            self.lines = fout.readlines()
        self.get_info()

    def get_info(self):
        """
        get the general information of opt run from opt run output file
        which is now stored in self.lines
        """
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
        for i in range(begin_final_coord_line+4, begin_final_coord_line+7):
            for j in range(3):
                self.cell.append(float(self.lines[i].split()[j]))
        for i in range(begin_final_coord_line+9, end_final_coord_line):
            self.atoms.append(Atom(self.lines[i].split()[0], float(self.lines[i].split()[1]), float(self.lines[i].split()[2]), float(self.lines[i].split()[3])))
   

    #
    
    def get_opt_params_and_run_info(self):
        """
        """
        self.run_info["iterations"] = []
        self.run_info["total-energies"] = []
        self.run_info["fermi-energies"] = []
        self.run_info["total-forces"] = []

        for line in self.lines:
            # if it is an empty line continue to next line
            if len(line.split()) == 0:
                continue
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
        cell = self.cell
        with open(xyz, 'w') as fout:
            fout.write("%d\n" % len(self.atoms))
            if self.run_type == "vc-relax":
                fout.write("cell: %f %f %f | %f %f %f | %f %f %f\n" % (cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7], cell[8]))
            else:
                fout.write("type of opt run: relax -> the cell is not changed, so go and find the original cell\n")
            for atom in self.atoms:
                fout.write("%s\t%f\t%f\t%f\n" % (atom.name, atom.x, atom.y, atom.z))
    
    def plot_run_info(self):
        """
        """
        plt.plot(self.run_info["iterations"])
        plt.title("Iterations per SCF")
        plt.xlabel("Scf scycles")
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

            fout.write("Total forces per SCF\n")
            fout.write("![Total forces per SCF](total-forces-per-scf.png)\n")


    def export(self):
        self.to_xyz()
        self.plot_run_info()
        self.markdown_report("OptimizationReport.md")