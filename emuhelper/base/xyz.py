#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import copy
import pymatgen as mg

from emuhelper.base.atom import Atom

"""
Usage:
"""


class base_xyz:
    """
    a representation of xyz structure
    usage:
        a = base_xyz()
        a.get_info(xyzfile)
    """
    def __init__(self):
        #self.file = xyz_f
        self.natom = 0
        self.nspecies = 0
        self.atoms = []
        self.specie_labels = dict()
        #self.get_info()
        #self.cell = self.get_cell()
        self.cell = None

    def get_info(self, xyz_f):
        """
        get info to construct the xyz from an xyz file
        """
        with open(xyz_f, 'r') as fin:
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
        self.cell = self.get_cell(xyz_f)

    def set_species_number(self):
        names = [self.atoms[x].name for x in range(self.natom)]
        species = set(names)
        species = list(species)
        species_with_order = {}
        for i in species:
            species_with_order[i] = mg.Element(i).number
        tmp = sorted(zip(species_with_order.values(), species_with_order.keys()))
        for i in range(len(tmp)):
            tmp[i] = list(tmp[i])
        for i in range(len(tmp)):
            tmp[i][0] = i + 1
        tmp = dict(tmp)
        self.specie_labels = dict(zip(tmp.values(), tmp.keys()))
        self.nspecies = len(self.specie_labels)


    def get_cell(self, xyz_f):
        """
        cell defined in xxx.xyz must be in format like this:
        cell: 4.08376 0.00000 0.00000 | 0.00000 4.00251 0.00000 | -0.05485 0.00000 8.16247
        """
        with open(xyz_f, 'r') as fin:
            fin.readline()
            line = fin.readline()
        return [float(line.split()[i]) for i in [1, 2, 3, 5, 6, 7, 9, 10, 11]]

    def update(self, newxyzfile):
        #self.file = newxyzfile
        self.natom = 0
        self.nspecies = 0
        self.atoms = []
        self.specie_labels = dict()
        #self.get_info()
        #self.cell = self.get_cell()
        self.get_info(newxyzfile)
    
    def build_supercell(self, n):
        """
        n: [n1, n2, n3]
        """
        self.origin_atoms = copy.deepcopy(self.atoms)
        self.origin_cell = copy.deepcopy(self.cell)
        self.origin_natom = copy.deepcopy(self.natom)
        #
        for i in range(3):
            self.cell[3*i+0] = n[i] * self.cell[3*i+0]
            self.cell[3*i+1] = n[i] * self.cell[3*i+1]
            self.cell[3*i+2] = n[i] * self.cell[3*i+2]
        #
        # clone the atoms to build supercell    
        for i in range(3):
            natom_now = len(self.atoms)
            for j in range(n[i] - 1):
                for atom in self.atoms[:natom_now]:
                    x = atom.x + float(j + 1) * self.origin_cell[3*i+0]
                    y = atom.y + float(j + 1) * self.origin_cell[3*i+1]
                    z = atom.z + float(j + 1) * self.origin_cell[3*i+2]
                    self.atoms.append(Atom(atom.name, x, y, z))
        self.natom = len(self.atoms)
    
    def to_xyz(self, fname):
        with open(fname, 'w') as fout:
            fout.write("%d\n" % self.natom)
            fout.write("cell: %f %f %f | %f %f %f | %f %f %f\n" % (self.cell[0], self.cell[1], self.cell[2], self.cell[3], self.cell[4], self.cell[5], self.cell[6], self.cell[7], self.cell[8]))
            for atom in self.atoms:
                fout.write("%s %f %f %f\n" % (atom.name, atom.x, atom.y, atom.z))



