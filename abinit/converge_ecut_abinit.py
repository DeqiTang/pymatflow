#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
    python converge_ecut.py xxx.xyz ecut_min ecut_max ecut_step
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

ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
ecut_max = int(sys.argv[3])
ecut_step = int(sys.argv[4])


xyz = XYZ()

base_project_name = "test"
if os.path.exists("./tmp-ecut"):
    shutil.rmtree("./tmp-ecut")
os.mkdir("./tmp-ecut")
os.chdir("./tmp-ecut")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
shutil.copyfile("../Li.psp8", "Li.psp8")
shutil.copyfile("../H.psp8", "H.psp8")

n_test = int((ecut_max - ecut_min) / ecut_step)
for i in range(n_test + 1):
    cutoff = int(ecut_min + i * ecut_step)
    rel_cutoff = cutoff / 3
    inp_name = "test-ecut-%d.in" % cutoff
    files_name = "test-ecut-%d.files" % cutoff
    with open(files_name, 'w') as fout:
        fout.write(inp_name)
        fout.write("\n")
        fout.write("test-ecut-%d.out\n" % cutoff)
        fout.write("test-ecut-%di\n" % cutoff)
        fout.write("test-ecut-%do\n" % cutoff)
        fout.write("temp\n")
        for element in xyz.specie_labels:
            fout.write("%s\n" % (element + ".psp8"))
        #
    with open(inp_name, 'w') as fout:
        fout.write("acell 10 10 10\n")
        fout.write("ntypat %d\n" % xyz.nspecies)
        fout.write("natom %d\n" % xyz.natom)

        fout.write("ecut %d\n" % cutoff)
        fout.write("kptopt 0\n")
        fout.write("nkpt 1\n")
        fout.write("nstep 100\n")
        fout.write("toldfe 1.0d-6\n")
        fout.write("diemac 2.0\n")
        fout.write("znucl ")
        for element in xyz.specie_labels:
            fout.write(str(mg.Element[element].number))
            fout.write(" ")
        fout.write("\n")
        fout.write("\n")
    xyz.to_abinit(inp_name)
# run the simulation
for i in range(n_test + 1):
    cutoff = int(ecut_min + i * ecut_step)
    files_name = "test-ecut-%d.files" % cutoff
    out_f_name = "test-ecut-%d.out" % cutoff
    os.system("abinit < %s > %s" % (files_name, out_f_name))

# analyse the result
for i in range(n_test + 1):
    cutoff = int(ecut_min + i * ecut_step)
    out_f_name = "test-ecut-%d.out" % cutoff
    os.system("cat %s | grep 'Etotal=' >> energy-ecut.data" % out_f_name)

ecut = [ ecut_min + i * ecut_step for i in range(n_test + 1)]
energy = []
with open("energy-ecut.data", 'r') as fin:
    for line in fin:
        energy.append(float(line.split()[2]))

import matplotlib.pyplot as plt

#for i in range(len(energy)):
#    energy[i] = energy[i] - 31
plt.plot(ecut, energy)
plt.show()
