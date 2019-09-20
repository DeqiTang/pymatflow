#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
    python converge_rel_cutoff.py xxx.xyz rel_cutoff_min rel_cutoff_max rel_cutoff_step cutoff
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

    def to_subsys(self, fname):
        cell = self.cell
        with open(fname, 'a') as fout:
            fout.write("\t&SUBSYS\n")
            for element in self.specie_labels:
                fout.write("\t\t&KIND %s\n" % element)
                fout.write("\t\t\tBASIS_SET SZV-GTH-PADE\n")
                fout.write("\t\t\tPOTENTIAL GTH-PADE\n")
                fout.write("\t\t&END KIND\n")
            fout.write("\t\t&CELL\n")
            #fout.write("\t\t\tABC %f %f %f\n" % (cell[0], cell[4], cell[8]))
            fout.write("\t\t\tA %f %f %f\n" % (cell[0], cell[1], cell[2]))
            fout.write("\t\t\tB %f %f %f\n" % (cell[3], cell[4], cell[5]))
            fout.write("\t\t\tC %f %f %f\n" % (cell[6], cell[7], cell[8]))
            fout.write("\t\t&END CELL\n")
            fout.write("\t\t&TOPOLOGY\n")
            fout.write("\t\t\tCOORD_FILE_FORMAT xyz\n")
            fout.write("\t\t\tCOORD_FILE_NAME %s\n" % sys.argv[1])
            fout.write("\t\t&END TOPOLOGY\n")
            fout.write("\t&END SUBSYS\n")
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

rel_cutoff_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
rel_cutoff_max = int(sys.argv[3])
rel_cutoff_step = int(sys.argv[4])
cutoff = int(sys.argv[5])

xyz = XYZ()

base_project_name = "test"

if os.path.exists("./tmp-rel-cutoff"):
    shutil.rmtree("./tmp-rel-cutoff")
os.mkdir("./tmp-rel-cutoff")
os.chdir("./tmp-rel-cutoff")
shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")

n_test = int((rel_cutoff_max - rel_cutoff_min) / rel_cutoff_step)
for i in range(n_test + 1):
    rel_cutoff = int(rel_cutoff_min + i * rel_cutoff_step)
    inp_name = "test-rel-cutoff-%d.inp" % rel_cutoff
    with open(inp_name, 'w') as fout:
        fout.write("&GLOBAL\n")
        fout.write("\tPROJECT\t%s\n" % (base_project_name + str(rel_cutoff)))
        fout.write("\tRUN_TYPE ENERGY_FORCE\n")
        fout.write("\tPRINT_LEVEL LOW\n")
        fout.write("&END GLOBAL\n")
        fout.write("\n")

        fout.write("&FORCE_EVAL\n")
        fout.write("\tMETHOD Quickstep\n")
    # subsys
    xyz.to_subsys(inp_name)
    # end subsys
    with open(inp_name, 'a') as fout:
        # dft
        fout.write("\t&DFT\n")
        fout.write("\t\tBASIS_SET_FILE_NAME BASIS_SET\n")
        fout.write("\t\tPOTENTIAL_FILE_NAME GTH_POTENTIALS\n")
        fout.write("\t\t&QS\n")
        fout.write("\t\t\tEPS_DEFAULT 1.0E-10\n")
        fout.write("\t\t&END QS\n")
        fout.write("\t\t&MGRID\n")
        fout.write("\t\t\tNGRIDS 4\n")
        fout.write("\t\t\tCUTOFF %d\n" % cutoff)
        fout.write("\t\t\tREL_CUTOFF %d\n" % rel_cutoff)
        fout.write("\t\t&END MGRID\n")
        fout.write("\t\t&XC\n")
        fout.write("\t\t\t&XC_FUNCTIONAL PADE\n")
        fout.write("\t\t\t&END XC_FUNCTIONAL\n")
        fout.write("\t\t&END XC\n")
        fout.write("\t\t&SCF\n")
        fout.write("\t\t\tSCF_GUESS ATOMIC\n")
        fout.write("\t\t\tEPS_SCF 1.0E-06\n")
        fout.write("\t\t\tMAX_SCF 300\n")
        fout.write("\t\t\t&DIAGONALIZATION ON\n")
        fout.write("\t\t\t\tALGORITHM STANDARD\n")
        fout.write("\t\t\t&END DIAGONALIZATION\n")
        fout.write("\t\t\t&MIXING T\n")
        fout.write("\t\t\t\tMETHOD BROYDEN_MIXING\n")
        fout.write("\t\t\t\tALPHA 0.4\n")
        fout.write("\t\t\t\tNBROYDEN 8\n")
        fout.write("\t\t\t&END MIXING\n")
        fout.write("\t\t&END SCF\n")
        fout.write("\t&END DFT\n")
        # end dft
        fout.write("&END FORCE_EVAL\n")

# run the simulation
for i in range(n_test + 1):
    rel_cutoff = int(rel_cutoff_min + i * rel_cutoff_step)
    inp_name = "test-rel-cutoff-%d.inp" % rel_cutoff
    out_f_name = "test-rel-cutoff-%d.out" % rel_cutoff
    os.system("cp2k.psmp -in %s > %s" % (inp_name, out_f_name))


# analyse the result
for i in range(n_test + 1):
    rel_cutoff = int(rel_cutoff_min + i * rel_cutoff_step)
    out_f_name = "test-rel-cutoff-%d.out" % rel_cutoff
    os.system("cat %s | grep 'Total energy:' >> energy-rel-cutoff.data" % out_f_name)

rel_cut = [rel_cutoff_min + i * rel_cutoff_step for i in range(n_test + 1)]
energy = []
with open("energy-rel-cutoff.data", 'r') as fin:
    for line in fin:
        energy.append(float(line.split()[2]))

import matplotlib.pyplot as plt

plt.plot(rel_cut, energy)
plt.show()
