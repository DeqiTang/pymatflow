#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
    python neb_abinit.py initial.xyz final.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential file
    for all the elements of the system is in the directory.
Note:
    参考: https://docs.abinit.org/topics/TransPath/
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
        self.cell = self.get_cell()

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
        cell = self.cell
        with open(fname, 'a') as fout:
            fout.write("acell 1 1 1\n") # scaling with 1 means no actually scaling of rprim by acell
            fout.write("rprim\n")
            fout.write("%f %f %f\n" % (cell[0], cell[1], cell[2]))
            fout.write("%f %f %f\n" % (cell[3], cell[4], cell[5]))
            fout.write("%f %f %f\n" % (cell[6], cell[7], cell[8]))

            fout.write("ntypat %d\n" % self.nspecies)
            fout.write("natom %d\n" % self.natom)
            fout.write("typat ")
            for atom in self.atoms:
                fout.write("%d " % self.specie_labels[atom.name])
            fout.write("\n")
            fout.write("znucl ")
            for element in self.specie_labels:
                fout.write(str(mg.Element[element].number))
                fout.write(" ")
            fout.write("\n")
            fout.write("\n")
            fout.write("xangst\n")
            for atom in self.atoms:
                fout.write("%f %f %f\n" % (atom.x, atom.y, atom.z))
            fout.write("\n")
    def to_xangst_lastimg(self, fname):
        with open(fname, 'a') as fout:
            fout.write("xangst_lastimg\n")
            for atom in self.atoms:
                fout.write("%f %f %f\n" % (atom.x, atom.y, atom.z))
            fout.write("\n")

    def get_cell(self):
        """
        cell defined in xxx.xyz must be in format like this:
        cell: 4.08376 0.00000 0.00000 | 0.00000 4.00251 0.00000 | -0.05485 0.00000 8.16247
        """
        with open(self.file, 'r') as fin:
            fin.readline()
            line = fin.readline()
        return [float(line.split()[i]) for i in [1, 2, 3, 5, 6, 7, 9, 10, 11]]

    def update(self, newxyzfile):
        self.file = newxyzfile
        self.natom = 0
        self.nspecies = 0
        self.atoms = []
        self.specie_labels = dict()
        self.get_info()
        self.cell = self.get_cell()

        

# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])
cutoff = 40

xyz_initial = XYZ(sys.argv[1])
xyz_final = XYZ(sys.argv[2])

base_project_name = "neb-transpath"
if os.path.exists("./tmp-neb"):
    shutil.rmtree("./tmp-neb")
os.mkdir("./tmp-neb")
os.chdir("./tmp-neb")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../Li.psp8", "Li.psp8")
#shutil.copyfile("../H.psp8", "H.psp8")
os.system("cp ../*.psp8 ./")


inp_name = "neb.in"
files_name = "neb.files"
with open(files_name, 'w') as fout:
    fout.write(inp_name)
    fout.write("\n")
    fout.write("neb.out\n")
    fout.write("nebi\n")
    fout.write("nebo\n")
    fout.write("temp\n")
    for element in xyz_initial.specie_labels:
        fout.write("%s\n" % (element + ".psp8"))
    #
with open(inp_name, 'w') as fout:
    fout.write("ecut %d\n" % cutoff)
    fout.write("kptopt 1\n")
    fout.write("ngkpt 1 1 1\n")
    fout.write("occopt 3\n") # fermi dirac smearing of occupation
    fout.write("nstep 100\n")
    fout.write("toldff 1.0d-6\n")
    fout.write("nband 10\n")
    fout.write("diemac 2.0\n")
    fout.write("nimage 12\n")
    fout.write("imgmov 5\n")
    fout.write("ntimimage 50\n")
    fout.write("tolimg 0.0001\n") # Tol. criterion (will stop when average energy of cells < tolimg)
    fout.write("dynimage 0 10*1 0\n")  # Keep first and last images fixed
    fout.write("fxcartfactor 1.0\n")
    fout.write("prtvolimg 0\n")
    fout.write("\n")
xyz_initial.to_abinit(inp_name)
xyz_final.to_xangst_lastimg(inp_name)


# run the simulation
#out_f_name = "geo-opt-calc.out.log"
os.system("abinit < %s" % (files_name))

# analyse the result

import matplotlib.pyplot as plt

