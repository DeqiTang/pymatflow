#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
    python phonon_abinit.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential file
    for all the elements of the system is in the directory.
Reference:
    https://docs.abinit.org/tutorial/rf2/
Note:
    目前此脚本进行的声子计算还未成功
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

base_project_name = "phonon-calc"
if os.path.exists("./tmp-phonon"):
    shutil.rmtree("./tmp-phonon")
os.mkdir("./tmp-phonon")
os.chdir("./tmp-phonon")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
shutil.copyfile("../Li.psp8", "Li.psp8")
shutil.copyfile("../H.psp8", "H.psp8")


# 1 Generation of a derivative database
inp_name_1 = "phonon-calc-1.in"
files_name_1 = "phonon-calc-1.files"
with open(files_name_1, 'w') as fout:
    fout.write(inp_name_1)
    fout.write("\n")
    fout.write("phonon-calc-1.out\n")
    fout.write("phonon-calc-1i\n")
    fout.write("phonon-calc-1o\n")
    fout.write("temp\n")
    for element in xyz.specie_labels:
        fout.write("%s\n" % (element + ".psp8"))
    #
with open(inp_name_1, 'w') as fout:
    # computation of the phonon spectrum
    fout.write("ndtset 10\n")
    # Set 1 : ground state self-consistency
    fout.write("getwfk1 0 # cancel default\n")
    fout.write("kptopt1 1 # automatic generation of k points, taking into account the symmetry\n")
    fout.write("nqpt1 0 # cancel default\n")
    fout.write("tolvrs1 1.0d-18 # SCF stopping criterion (modify default)\n")
    fout.write("rfphon1 0 # cancel default\n")
    fout.write("nqpt 1\n") #  Q vectors for all datasets
    fout.write("qpt2   0.00000000E+00  0.00000000E+00  0.00000000E+00\n")
    fout.write("qpt3   0.00000000E+00  0.00000000E+00  0.00000000E+00\n")
    fout.write("qpt4   2.50000000E-01  0.00000000E+00  0.00000000E+00\n")
    fout.write("qpt5   5.00000000E-01  0.00000000E+00  0.00000000E+00\n")
    fout.write("qpt6   2.50000000E-01  2.50000000E-01  0.00000000E+00\n")
    fout.write("qpt7   5.00000000E-01  2.50000000E-01  0.00000000E+00\n")
    fout.write("qpt8  -2.50000000E-01  2.50000000E-01  0.00000000E+00\n")
    fout.write("qpt9   5.00000000E-01  5.00000000E-01  0.00000000E+00\n")
    fout.write("qpt10 -2.50000000E-01  5.00000000E-01  2.50000000E-01\n")
    fout.write("\n")
    # Set 2 : Response function calculation of d/dk wave function 
    fout.write("iscf2 -3\n")
    fout.write("kptopt2 2\n")
    fout.write("rfphon2 0 # cancel default\n")
    fout.write("rfelfd2 2 \n")
    fout.write("tolwfr2 1.0d-22\n")
    fout.write("\n")
    # Set 3 : Response function calculation of Q=0 phonons and electric field pert.
    fout.write("getddk3 2\n")
    fout.write("kptopt3 2\n")
    fout.write("rfelfd3 3\n")
    fout.write("\n")
    # Sets 4-10 : Finite-wave-vector phonon calculations (defaults for all datasets)
    fout.write("getwfk 1\n")
    fout.write("kptopt 3\n")
    fout.write("rfphon 1 # Do phonon response\n")
    fout.write("rfatpol 1 2\n")
    fout.write("rfdir 1 1 1\n")
    fout.write("tolvrs 1.0d-8\n")
    
    # definition of general setup
    fout.write("acell 4.000175 4.000175 4.000175\n")
    fout.write("scalecart 3*1\n") # 一定要注意scalecart的重要性, 设置大了可能会耗尽内存, 但对于周期体系又和重要
    fout.write("ntypat %d\n" % xyz.nspecies)
    fout.write("natom %d\n" % xyz.natom)
    fout.write("ecut %d\n" % cutoff)
    fout.write("ngkpt 4 4 4\n") # 似乎k点需要设置为偶数, 设置为奇数比如 3 3 3 计算会报错,说k点不对称
    fout.write("nstep 50\n")
    fout.write("diemac 2.0\n")
    fout.write("znucl ")
    for element in xyz.specie_labels:
        fout.write(str(mg.Element[element].number))
        fout.write(" ")
    fout.write("\n")
    fout.write("\n")
xyz.to_abinit(inp_name_1)

# run the simulation
out_f_name_1 = "phonon-calc-1.out"
os.system("abinit < %s > %s" % (files_name_1, out_f_name_1))

# 2 Manipulation of the derivative databases (the MRGDDB utility)
inp_name_3 = "phonon-calc-3.in"
with open(inp_name_3, 'w') as fout:
    fout.write("phonon-calc-3.ddb.out\n")
    fout.write("System phonons on 3 3 3 mesh\n")
    fout.write("8\n")
    fout.write("phonon-calc-1o_DS3_DDB\n")
    fout.write("phonon-calc-1o_DS4_DDB\n")
    fout.write("phonon-calc-1o_DS5_DDB\n")
    fout.write("phonon-calc-1o_DS6_DDB\n")
    fout.write("phonon-calc-1o_DS7_DDB\n")
    fout.write("phonon-calc-1o_DS8_DDB\n")
    fout.write("phonon-calc-1o_DS9_DDB\n")
    fout.write("phonon-calc-1o_DS10_DDB\n")
os.system("mrgddb < %s" % inp_name_3)


# 3 Analysis of the derivative databases


# 4 The computation of interatomic force constants
inp_name_4 = "phonon-calc-4.in"
files_name_4 = "phonon-calc-4.files"
with open(files_name_4, 'w') as fout:
    fout.write("phonon-calc-4.in\n")
    fout.write("phonon-calc-4.out\n")
    fout.write("phonon-calc-3.ddb.out\n")
    fout.write("phonon-calc-4_band2eps\n")
    fout.write("phonon-calc_dummy1\n")
    fout.write("phonon-calc_dummy2\n")
    fout.write("phonon-calc_dummy3\n")
with open(inp_name_4, 'w') as fout:
    fout.write("ifcflag 1\n")
    fout.write("brav 1\n") # Bravais Lattice : 1-S.C., 2-F.C., 3-B.C., 4-Hex.)
    fout.write("ngqpt 4 4 4\n")
    fout.write("nqshft 1\n")
    fout.write("q1shft 3*0.0\n")
    fout.write("chneut 1\n")
    # Interatomic force constant info
    fout.write("dipdip 1\n")
    fout.write("ifcana 1\n")
    fout.write("ifcout 20\n")
    fout.write("natifc 1\n")
    fout.write("atifc 1\n")
    fout.write("\n")
os.system("anaddb < %s > phonon-calc-4.log" % (files_name_4))

"""
# 5 Computation of phonon band structures with efficient interpolation
inp_name_5 = "phonon-calc-5.in"
files_name_5 = "phonon-calc-5.files"
with open(files_name_5, 'w') as fout:
    fout.write("phonon-calc-5.in\n")
    fout.write("phonon-calc-5.out\n")
    fout.write("phonon-calc-3.ddb.out\n")
    fout.write("phonon-calc-5_band2eps\n")
    fout.write("phonon-calc_dummy1\n")
    fout.write("phonon-calc_dummy2\n")
    fout.write("phonon-calc_dummy3\n")
with open(inp_name_5, 'w') as fout:
    fout.write("ifcflag 1\n")
    fout.write("ifcout 0\n")
    fout.write("brav 1\n") # Bravais Lattice : 1-S.C., 2-F.C., 3-B.C., 4-Hex
    fout.write("ngqpt 4 4 4\n") # Monkhorst-Pack indices
    fout.write("nqshft 1\n")
    fout.write("q1shft 3*0.0\n")
    fout.write("chneut 1\n")
    fout.write("dipdip 1\n")
    fout.write("eivec 4\n")
    fout.write("nph1l 71\n") # number of phonons in list 1   
    fout.write("qph1l   0.0000  0.0000  0.0000   1.0\n")    # !(gamma point)
    fout.write("0.0375  0.0375  0.0750   1.0\n")
    fout.write("0.0750  0.0750  0.1500   1.0\n")
    fout.write("0.1125  0.1125  0.2250   1.0\n")
    fout.write("0.1500  0.1500  0.3000   1.0\n")
    fout.write("0.1875  0.1875  0.3750   1.0\n")
    fout.write("0.2250  0.2250  0.4500   1.0\n")
    fout.write("0.2625  0.2625  0.5250   1.0\n")
    fout.write("0.3000  0.3000  0.6000   1.0\n")
    fout.write("0.3375  0.3375  0.6750   1.0\n")
    fout.write("0.3750  0.3750  0.7500   1.0\n")   # !(K point)
    fout.write("0.3875  0.3875  0.7750   1.0\n")
    fout.write("0.4000  0.4000  0.8000   1.0\n")
    fout.write("0.4125  0.4125  0.8250   1.0\n")
    fout.write("0.4250  0.4250  0.8500   1.0\n")
    fout.write("0.4375  0.4375  0.8750   1.0\n")
    fout.write("0.4500  0.4500  0.9000   1.0\n")
    fout.write("0.4625  0.4625  0.9250   1.0\n")
    fout.write("0.4750  0.4750  0.9500   1.0\n")
    fout.write("0.4875  0.4875  0.9750   1.0\n")
    fout.write("0.5000  0.5000  1.0000   1.0\n")    # !(X point)
    fout.write("0.5500  0.5500  1.0000   1.0\n")
    fout.write("0.6000  0.6000  1.0000   1.0\n")
    fout.write("0.6500  0.6500  1.0000   1.0\n")
    fout.write("0.7000  0.7000  1.0000   1.0\n")
    fout.write("0.7500  0.7500  1.0000   1.0\n")
    fout.write("0.8000  0.8000  1.0000   1.0\n")
    fout.write("0.8500  0.8500  1.0000   1.0\n")
    fout.write("0.9000  0.9000  1.0000   1.0\n")
    fout.write("0.9500  0.9500  1.0000   1.0\n")
    fout.write("1.0000  1.0000  1.0000   1.0\n")    # !(gamma point)
    fout.write("0.9500  0.9500  0.9500   1.0\n")
    fout.write("0.9000  0.9000  0.9000   1.0\n")
    fout.write("0.8500  0.8500  0.8500   1.0\n")
    fout.write("0.8000  0.8000  0.8000   1.0\n")
    fout.write("0.7500  0.7500  0.7500   1.0\n")
    fout.write("0.7000  0.7000  0.7000   1.0\n")
    fout.write("0.6500  0.6500  0.6500   1.0\n")
    fout.write("0.6000  0.6000  0.6000   1.0\n")
    fout.write("0.5500  0.5500  0.5500   1.0\n")
    fout.write("0.5000  0.5000  0.5000   1.0\n")     # !(L point)
    fout.write("0.5000  0.4500  0.5000   1.0\n")
    fout.write("0.5000  0.4000  0.5000   1.0\n")
    fout.write("0.5000  0.3500  0.5000   1.0\n")
    fout.write("0.5000  0.3000  0.5000   1.0\n")
    fout.write("0.5000  0.2500  0.5000   1.0\n")
    fout.write("0.5000  0.2000  0.5000   1.0\n")
    fout.write("0.5000  0.1500  0.5000   1.0\n")
    fout.write("0.5000  0.1000  0.5000   1.0\n")
    fout.write("0.5000  0.0500  0.5000   1.0\n")
    fout.write("0.5000  0.0000  0.5000   1.0\n")     # !(X point)
    fout.write("0.5000  0.0250  0.5250   1.0\n")
    fout.write("0.5000  0.0500  0.5500   1.0\n")
    fout.write("0.5000  0.0750  0.5750   1.0\n")
    fout.write("0.5000  0.1000  0.6000   1.0\n")
    fout.write("0.5000  0.1250  0.6250   1.0\n")
    fout.write("0.5000  0.1500  0.6500   1.0\n")
    fout.write("0.5000  0.1750  0.6750   1.0\n")
    fout.write("0.5000  0.2000  0.7000   1.0\n")
    fout.write("0.5000  0.2250  0.7250   1.0\n")
    fout.write("0.5000  0.2500  0.7500   1.0\n")     # !(W point)
    fout.write("0.5000  0.2750  0.7250   1.0\n")
    fout.write("0.5000  0.3000  0.7000   1.0\n")
    fout.write("0.5000  0.3250  0.6750   1.0\n")
    fout.write("0.5000  0.3500  0.6500   1.0\n")
    fout.write("0.5000  0.3750  0.6250   1.0\n")
    fout.write("0.5000  0.4000  0.6000   1.0\n")
    fout.write("0.5000  0.4250  0.5750   1.0\n")
    fout.write("0.5000  0.4500  0.5500   1.0\n")
    fout.write("0.5000  0.4750  0.5250   1.0\n")
    fout.write("0.5000  0.5000  0.5000   1.0\n")     # !(L point)
    fout.write("nph2l 1\n") # number of directions in list 2
    fout.write("qph2l 1.0 0.0 0.0 0.0  0.0\n")

os.system("anaddb < %s >phonon-calc-5.log" % (files_name_5))


files_name_6 = "phonon-calc-6.in"
inp_name_6 = "phonon-calc-6.in"
with open(files_name_6, 'w') as fout:
    fout.write("phonon-calc-6.in\n")
    fout.write("phonon-calc-6.out.eps\n")
    fout.write("phonon-calc-6.out_B2EPS.freq\n")
    fout.write("no\n")

with open(inp_name_6, 'w') as fout:
    fout.write("natom 2\n")
    fout.write("min 0.0 max 400.0 ngrad 8\n")
    fout.write("cunit 1\n")
    fout.write("nlines 7\n")
    fout.write("qpoint_name gamma K X gamma L X W L\n")
    fout.write("nqline 10 10 10 10 10 10 11\n")
    fout.write("scale 1.06066017  0.35355339  1.0  0.86602540  0.86602540  0.5  0.70710678\n")
    fout.write("red 0 0\n")
    fout.write("green 0 0\n")
    fout.write("blue 0 0\n")
#
os.system("band2eps < %s > phponon-calc-6.log" % (files_name_6))
# if you don't have gv: trizen -S --noconfirm gv
os.system("gv phonon-calc-6.out.eps")
os.system("abiopen.py phonon-calc-5.out_PHBST.nc --expose --seaborn=talk\n")

os.system("cp phonon-calc-3.ddb.out phonon-calc-3_DDB")
os.system("abiview.py ddb phonon-calc-3_DDB -sns=talk")

# 6 Thermodynamical properties
# analyse the result
"""
