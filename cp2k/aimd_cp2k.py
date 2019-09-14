#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
    python aimd_cp2k.py xxx.xyz 
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

    def to_subsys(self, fname):
        with open(fname, 'a') as fout:
            fout.write("\t&SUBSYS\n")
            for element in self.specie_labels:
                fout.write("\t\t&KIND %s\n" % element)
                fout.write("\t\t\tBASIS_SET DZV-GTH-PADE\n")
                fout.write("\t\t\tPOTENTIAL GTH-PADE\n")
                fout.write("\t\t&END KIND\n")
            fout.write("\t\t&CELL\n")
            fout.write("\t\t\tABC 10 10 10\n")
            fout.write("\t\t&END CELL\n")
            fout.write("\t\t&TOPOLOGY\n")
            fout.write("\t\t\tCOORD_FILE_FORMAT xyz\n")
            fout.write("\t\t\tCOORD_FILE_NAME %s\n" % sys.argv[1])
            fout.write("\t\t&END TOPOLOGY\n")
            fout.write("\t&END SUBSYS\n")
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

#cutoff_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#cutoff_max = int(sys.argv[3])
#cutoff_step = int(sys.argv[4])
#rel_cutoff = int(sys.argv[5])
cutoff = 60
rel_cutoff = 30

xyz = XYZ()

base_project_name = "aimd"

if os.path.exists("./tmp-aimd"):
    shutil.rmtree("./tmp-aimd")
os.mkdir("./tmp-aimd")
os.chdir("./tmp-aimd")
shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")

inp_name = "aimd.inp"
with open(inp_name, 'w') as fout:
    fout.write("&GLOBAL\n")
    fout.write("\tPROJECT\t%s\n" % (base_project_name ))
    fout.write("\tRUN_TYPE MD\n")
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
    fout.write("\t\t\tMAX_SCF 200\n")
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

    # MOTION
    fout.write("&MOTION\n")
    fout.write("\t&MD\n")
    fout.write("\t\tENSEMBLE NVT\n")
    fout.write("\t\tSTEPS 100\n")
    fout.write("\t\tTIMESTEP 0.5\n")
    fout.write("\t\t&THERMOSTAT\n")
    fout.write("\t\t\tTYPE NOSE\n")
    fout.write("\t\t\t&NOSE\n")
    fout.write("\t\t\t\tTIMECON 10.0\n")
    fout.write("\t\t\t&END NOSE\n")
    fout.write("\t\t&END THERMOSTAT\n")
    fout.write("\t\tTEMPERATURE 300.0\n")
    fout.write("\t&END MD\n")
    fout.write("\t&PRINT\n")
    fout.write("\t\t&RESTART\n")
    fout.write("\t\t\t&EACH\n")
    fout.write("\t\t\t\tMD 0\n")
    fout.write("\t\t\t&END EACH\n")
    fout.write("\t\t&END RESTART\n")
    fout.write("\t&END PRINT\n")
    fout.write("&END MOTION\n")
    fout.write("\n")

# run the simulation
out_f_name = "aimd.out"
os.system("cp2k.psmp -in %s > %s" % (inp_name, out_f_name))


# analyse the result

import matplotlib.pyplot as plt
