#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
    python bands_abinit.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential file
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

    def to_abinit(self, fname):
        with open(fname, 'a') as fout:
            fout.write("typat ")
            for atom in self.atoms:
                fout.write("%d " % self.specie_labels[atom.name])
            fout.write("\n")
            fout.write("xangst\n")
            for atom in self.atoms:
                fout.write("%d %d %d\n" % (atom.x, atom.y, atom.z))
            fout.write("\n")

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

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])
cutoff = 100

xyz = XYZ()

base_project_name = "bands-calc"
if os.path.exists("./tmp-bands"):
    shutil.rmtree("./tmp-bands")
os.mkdir("./tmp-bands")
os.chdir("./tmp-bands")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
shutil.copyfile("../Li.psp8", "Li.psp8")
shutil.copyfile("../H.psp8", "H.psp8")

inp_name = "bands-calc.in"
files_name = "bands-calc.files"
with open(files_name, 'w') as fout:
    fout.write(inp_name)
    fout.write("\n")
    fout.write("bands-calc.out\n")
    fout.write("bands-calci\n")
    fout.write("bands-calco\n")
    fout.write("temp\n")
    for element in xyz.specie_labels:
        fout.write("%s\n" % (element + ".psp8"))
    #
with open(inp_name, 'w') as fout:
    fout.write("ndtset 2\n")
    # dataset 1: scf calculation
    fout.write("# Dataset 1: scf\n")
    fout.write("kptopt1 1\n")
    fout.write("nshiftk1 4\n")
    fout.write("shiftk1 0.5 0.5 0.5\n")
    fout.write("  0.5 0.0 0.0\n")
    fout.write("  0.0 0.5 0.0\n")
    fout.write("  0.0 0.0 0.5\n")
    fout.write("ngkpt1 4 4 4\n")
    fout.write("prtden1 1\n")
    fout.write("toldfe1 1.0d-6\n")
    # dataset 2: band calculation
    fout.write("iscf2 -2\n")
    fout.write("getden2 -1\n")
    fout.write("kptopt2 -3\n")
    fout.write("nband2 8\n")
    fout.write("ndivk2 10 12 17\n")
    fout.write("kptbounds2 0.5 0.0 0.0 # L point\n")
    fout.write("0.0 0.0 0.0 # Gamma\n")
    fout.write("0.0 0.5 0.5 # X\n")
    fout.write("1.0 1.0 1.0 # Gamma in another cell\n")
    fout.write("tolwfr2 1.0d-12\n")
    fout.write("enunit2 1\n")
    fout.write("\n")
    # define unit cell
    fout.write("acell 4.000175 4.000175 4.000175\n")
    fout.write("ntypat %d\n" % xyz.nspecies)
    fout.write("natom %d\n" % xyz.natom)

    fout.write("ecut %d\n" % cutoff)
    fout.write("occopt 3\n")
    #fout.write("prtdos 1\n") # prtdos = 1, 2, 3
    fout.write("nstep 100\n")
    fout.write("diemac 2.0\n")
    fout.write("znucl ")
    for element in xyz.specie_labels:
        fout.write(str(mg.Element[element].number))
        fout.write(" ")
    fout.write("\n")
    fout.write("\n")
xyz.to_abinit(inp_name)

# run the simulation
out_f_name = "bands-calc.out"
os.system("abinit < %s > %s" % (files_name, out_f_name))

# analyse the result

import matplotlib.pyplot as plt

