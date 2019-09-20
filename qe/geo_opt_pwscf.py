#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import re
import pymatgen as mg


"""
Usage:
    python geo_opt_pwscf.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential
    file for all the elements in the system is in the directory.
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

    def to_pwscf(self, fname):
        cell = self.cell
        with open(fname, 'a') as fout:
            fout.write("ATOMIC_SPECIES\n")
            for element in self.specie_labels:
                tmp = os.listdir("../")
                pseudo_file = ""
                for f in tmp:
                    match = re.match(element, f)
                    if match is not None:
                        pseudo_file = match.string
                        break
                fout.write("%s %f %s\n" % (element, mg.Element(element).atomic_mass, pseudo_file))

            fout.write("CELL_PARAMETERS (angstrom)\n")
            fout.write("%f %f %f\n" % (cell[0], cell[1], cell[2]))
            fout.write("%f %f %f\n" % (cell[3], cell[4], cell[5]))
            fout.write("%f %f %f\n" % (cell[6], cell[7], cell[8]))
            fout.write("\n")
            fout.write("ATOMIC_POSITIONS (crystal)\n")
            for atom in self.atoms:
                fout.write("%s\t%f\t%f\t%f\n" % (atom.name, atom.x, atom.y, atom.z))
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
        self.set_species_number()


        

# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])

ecut_wfc = 80
ecut_rho = ecut_wfc * 2

xyz = XYZ()

title = "Geometric Optimization"

base_prefix = "knn"

if os.path.exists("./tmp-geo-opt"):
    shutil.rmtree("./tmp-geo-opt")
os.mkdir("./tmp-geo-opt")
os.chdir("./tmp-geo-opt")

inp_name = "geo-opt.in"
with open(inp_name, 'w') as fout:
    fout.write("&control\n")
    fout.write("calculation = 'scf'\n")
    fout.write("title = '%s'\n" % title)
    fout.write("prefix = '%s'\n" % (base_prefix))
    fout.write("restart_mode = 'from_scratch'\n")
    fout.write("nstep = 300\n")
    fout.write("outdir = '%s'\n" % ("./tmp"))
    fout.write("pseudo_dir = '../'\n")
    fout.write("wf_collect = .true.\n")
    fout.write("tstress = .true.\n")
    fout.write("tprnfor = .true.\n")
    fout.write("/\n")
    fout.write("\n")

    fout.write("&system\n")
    fout.write("ibrav = 0\n")
    fout.write("nat = %d\n" % xyz.natom)
    fout.write("ntyp = %d\n" % xyz.nspecies)
    fout.write("nspin = 1\n")
    #fout.write("nspin = 2\n")
    #fout.write("starting_magnetization(1) = 1\n")
    #fout.write("starting_magnetization(2) = 1\n")
    fout.write("ecutwfc = %d\n" % ecut_wfc)
    fout.write("ecutrho = %d\n" % ecut_rho)
    fout.write("input_DFT = 'PBE'\n")
    fout.write("occupations = 'smearing'\n")
    fout.write("degauss = 1.0d-4\n")
    fout.write("smearing = 'marzari-vanderbilt'\n")
    fout.write("/\n")
    fout.write("\n")

    fout.write("&electrons\n")
    fout.write("electron_maxstep = 300\n")
    #fout.write("conv_thr = 1.0d-10\n")
    fout.write("conv_thr = 1.0d-5\n")
    fout.write("mixing_mode = 'plain'\n")
    fout.write("mixing_beta = 0.3d0\n")
    fout.write("scf_must_converge = .true.\n")
    fout.write("/\n")
    fout.write("\n")

    fout.write("K_POINTS automatic\n")
    fout.write("2 2 2 0 0 0\n")
    fout.write("\n")
xyz.to_pwscf(inp_name)


# run the simulation
inp_name = "geo-opt.in"
out_f_name = "geo-opt.out"
os.system("pw.x < %s > %s" % (inp_name, out_f_name))


# analyse the result

import matplotlib.pyplot as plt

