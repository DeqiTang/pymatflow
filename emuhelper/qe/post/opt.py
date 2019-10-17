#!/usr/bin/env python
# _*_ coding: utf-8 _*_

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
        
        self.get_opt_params()

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
    #
    
    def get_opt_params(self):
        """
        """
        for line in self.lines:
            if line.split()[0] == "kinetic-energy":
                self.opt_params["ecutwfc"] = int(float(line.split()[3]))
            if line.split()[0] == "convergence threshold":
                self.opt_params["conv_thr"] = float(line.split()[3])
            if line.split()[0] == "mixing" and line.split()[1] == 'beta':
                self.opt_params["mixing_beta"] = float(line.split()[3])
            if line.split()[0] == "number" and line.split()[2] == 'k':
                self.opt_params["degauss"] = float(line.split()[9])
