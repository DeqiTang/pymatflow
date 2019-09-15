#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
    python single_point_dalton.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.
"""


class Atom:
    """
    a representation of atom with xyz coordinates
    """
    def __init__(self, name=None, x=0, y=0, z=0):
        self.name = name
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
    def set_name(self, name):
        self.name = name
    def set_x(self, x):
        self.x = float(x)
    def set_y(self, y):
        self.y = float(y)
    def set_z(self, z):
        self.z = float(z)


class XYZ:
    """
    a representation of xyz file
    """
    def __init__(self, xyz_f=sys.argv[1]):
        self.file = xyz_f
        self.natom = 0
        self.nspecies = 0
        self.atoms = []
        self.specie_labels = dict()
        self.get_info()

    def get_info(self):
        with open(self.file, 'r') as fin:
            self.natom = int(fin.readline())
            fin.readline()
            i = 0
            while i < self.natom:
                line = fin.readline()
                atom = Atom(line.split()[0], float(line.split()[1]), float(line.split()[2]), float(line.split()[3]))
                self.atoms.append(atom)
                i += 1
        self.set_species_number()

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

    def to_dalton(self, fname):
        with open(fname, 'a') as fout:
            fout.write("Atomtypes=%d\n" % self.nspecies)
            for element in self.specie_labels:
                charge = float(mg.Element[element].number)
                num_atoms = 0
                for atom in self.atoms:
                    if atom.name == element:
                        num_atoms += 1
                fout.write("Charge=%f Atoms=%d\n" % (charge, num_atoms))
                for atom in self.atoms:
                    if atom.name == element:
                        fout.write("%s %f %f %f\n" % (atom.name, atom.x, atom.y, atom.z))
            # end 

    def update(self, newxyzfile):
        self.file = newxyzfile
        self.natom = 0
        self.nspecies = 0
        self.atoms = []
        self.specie_labels = dict()
        self.get_info()
        self.set_species_number()


        

# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]


xyz = XYZ()

#base_project_name = "test"

if os.path.exists("./tmp-single-point"):
    shutil.rmtree("./tmp-single-point")
os.mkdir("./tmp-single-point")
os.chdir("./tmp-single-point")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])


comment_1 = "xxx"
comment_2 = "xxx"

mol_name = "single-point.mol"
dal_name = "single-point.dal"
with open(mol_name, 'w') as fout:
    fout.write("BASIS\n")
    fout.write("6-31G**\n")
    fout.write("%s\n" % comment_1)
    fout.write("%s\n" % comment_2)
xyz.to_dalton(mol_name)

with open(dal_name, 'w') as fout:
    fout.write("**DALTON INPUT\n")
    fout.write(".RUN PROPERTIES\n")
    fout.write("**WAVE FUNCTIONS\n")
    fout.write(".HF\n")
    fout.write(".MP2\n")
    fout.write("*SCF INPUT\n")
    fout.write(".DOUBLY OCCUPIED\n")
    fout.write(" 3 1 1 0\n")
    fout.write("*CONFIGURATION INPUT\n")
    fout.write("**START\n")
    fout.write("**PROPERTIES\n")
    fout.write(".DIPGRA\n")
    fout.write(".VIBANA\n")
    fout.write("**END OF DALTON INPUT\n")

# run the simulation
#os.system("dalton -mol single-point -dal single-point")
os.system("dalton -N 2 -mol %s -dal %s" % (mol_name, dal_name))
