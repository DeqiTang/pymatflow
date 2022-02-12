"""
base xyz structure class, a representation of xyz structure
"""
import numpy as np
import sys
import os
import shutil
import copy
from pymatflow.base.element import element
from pymatflow.base.atom import Atom

"""
Usage:
"""


class BaseXyz:
    """
    a representation of xyz structure
    usage:
        a = base_xyz()
        a.get_xyz(xyzfile)
    """
    def __init__(self):
        self.file = None
        self.natom = 0
        self.nspecies = 0
        self.atoms = []
        self.cell = None
        self.specie_labels = dict()

    def get_xyz(self, xyzfile):
        """
            get information to construct the structure from an xyz file
        """
        #self.file = xyzfile
        self.file = os.path.abspath(xyzfile)
        # now self.file is an absolute path, this is much easier for later usage
        with open(self.file, 'r') as fin:
            self.natom = int(fin.readline())
            fin.readline()
            i = 0
            while i < self.natom:
                line = fin.readline()
                atom = Atom(line.split()[0], float(line.split()[1]), float(line.split()[2]), float(line.split()[3]))
                # get information about the fix of this atom for opt and md
                if len(line.split()) == 4:
                    atom.fix = [False, False, False]
                elif line.split()[4][0] == '#':
                    atom.fix = [False, False, False] # the first char after coordinates is # , so cannpt set T or F
                else:
                    for j in range(3):
                        if line.split()[j+4] == 'T':
                            atom.fix[j] = True
                        elif line.split()[j+4] == 'F':
                            atom.fix[j] = False
                        else:
                            print("===============================\n")
                            print("warning: while read xyz file\n")
                            print("can only set T or F after coords\n")
                            sys.exit(1)
                # end get the information about the fix of this atom for opt and md
                self.atoms.append(atom)
                i += 1
        self.set_species_number()
        self.cell = self.get_cell(self.file)

    def set_species_number(self):
        names = [self.atoms[x].name for x in range(self.natom)]
        species = set(names)
        species = list(species)
        species_with_order = {}
        for i in species:
            species_with_order[i] = element[i].number
        tmp = sorted(zip(species_with_order.values(), species_with_order.keys()))
        for i in range(len(tmp)):
            tmp[i] = list(tmp[i])
        for i in range(len(tmp)):
            tmp[i][0] = i + 1
        tmp = dict(tmp)
        self.specie_labels = dict(zip(tmp.values(), tmp.keys()))
        self.nspecies = len(self.specie_labels)


    def get_cell(self, xyzfile):
        """
            cell defined in xxx.xyz must be in format like this:
            cell: 4.08376 0.00000 0.00000 | 0.00000 4.00251 0.00000 | -0.05485 0.00000 8.16247
        """
        with open(xyzfile, 'r') as fin:
            fin.readline()
            line = fin.readline()
        #return [float(line.split()[i]) for i in [1, 2, 3, 5, 6, 7, 9, 10, 11]]
        cell = []
        cell.append([float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])
        cell.append([float(line.split()[5]), float(line.split()[6]), float(line.split()[7])])
        cell.append([float(line.split()[9]), float(line.split()[10]), float(line.split()[11])])

        return cell

    def update(self, xyzfile):
        self.file = xyzfile
        self.natom = 0
        self.nspecies = 0
        self.atoms = []
        self.specie_labels = dict()
        self.get_xyz(xyzfile=xyzfile)

    def to_xyz(self, fname):
        with open(fname, 'w') as fout:
            fout.write("%d\n" % len(self.atoms))
            #fout.write("cell: %f %f %f | %f %f %f | %f %f %f\n" % (self.cell[0], self.cell[1], self.cell[2], self.cell[3], self.cell[4], self.cell[5], self.cell[6], self.cell[7], self.cell[8]))
            fout.write("cell: %f %f %f | %f %f %f | %f %f %f\n" % (self.cell[0][0], self.cell[0][1], self.cell[0][2], self.cell[1][0], self.cell[1][1], self.cell[1][2], self.cell[2][0], self.cell[2][1], self.cell[2][2]))
            for atom in self.atoms:
                fout.write("%s %f %f %f\n" % (atom.name, atom.x, atom.y, atom.z))

    def to_xyz_file(self, fname):
        with open(fname, 'w') as fout:
            fout.write("%d\n" % len(self.atoms))
            #fout.write("cell: %f %f %f | %f %f %f | %f %f %f\n" % (self.cell[0], self.cell[1], self.cell[2], self.cell[3], self.cell[4], self.cell[5], self.cell[6], self.cell[7], self.cell[8]))
            fout.write("cell: %f %f %f | %f %f %f | %f %f %f\n" % (self.cell[0][0], self.cell[0][1], self.cell[0][2], self.cell[1][0], self.cell[1][1], self.cell[1][2], self.cell[2][0], self.cell[2][1], self.cell[2][2]))
            for atom in self.atoms:
                fout.write("%s %f %f %f\n" % (atom.name, atom.x, atom.y, atom.z))
