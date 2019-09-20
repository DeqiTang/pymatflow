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
    python converge_ecutwfc_pwscf.py xxx.xyz ecut_wfc_min ecut_wfc_max ecut_wfc_step
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential
    file for all the elements in the system is in the directory.
Note:
    here while finding the best ecutwfc I always set the the ecutrho to the default
    value: 4 * ecutwfc
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

ecut_wfc_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
ecut_wfc_max = int(sys.argv[3])
ecut_wfc_step = int(sys.argv[4])


xyz = XYZ()

title = "Test ecutwfc"
base_prefix = "knn"

if os.path.exists("./tmp-ecutwfc"):
    shutil.rmtree("./tmp-ecutwfc")
os.mkdir("./tmp-ecutwfc")
os.chdir("./tmp-ecutwfc")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")

n_test = int((ecut_wfc_max - ecut_wfc_min) / ecut_wfc_step)
for i in range(n_test + 1):
    ecut_wfc = int(ecut_wfc_min + i * ecut_wfc_step)
    ecut_rho = ecut_wfc * 4 # using default value for ecut_rho: 4 * ecutwfc
    inp_name = "test-ecutwfc-%d.in" % ecut_wfc
    with open(inp_name, 'w') as fout:
        fout.write("&control\n")
        fout.write("calculation = 'scf'\n")
        fout.write("title = '%s'\n" % title)
        fout.write("prefix = '%s'\n" % (base_prefix + str(ecut_wfc)))
        fout.write("restart_mode = 'from_scratch'\n")
        fout.write("nstep = 300\n")
        fout.write("outdir = '%s'\n" % ("./tmp-" + str(ecut_wfc)))
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
        fout.write("ecutrho = %d \n" % ecut_rho)  # default value: 4 * ecutwfc
        fout.write("input_DFT = 'PBE'\n")
        fout.write("occupations = 'fixed'\n") # smearing, tetrahedra, fixed
        #fout.write("degauss = 1.0d-4\n") # default: 0
        #fout.write("smearing = 'gaussian'\n") # default is gaussian, and you can use fermi-dirac
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
        fout.write("1 1 1 0 0 0\n")
        fout.write("\n")
    xyz.to_pwscf(inp_name)


# run the simulation
for i in range(n_test + 1):
    ecut_wfc = int(ecut_wfc_min + i * ecut_wfc_step)
    inp_name = "test-ecutwfc-%d.in" % ecut_wfc
    out_f_name = "test-ecutwfc-%d.out" % ecut_wfc
    os.system("pw.x < %s > %s" % (inp_name, out_f_name))


# analyse the result
for i in range(n_test + 1):
    ecut_wfc = int(ecut_wfc_min + i * ecut_wfc_step)
    out_f_name = "test-ecutwfc-%d.out" % ecut_wfc
    os.system("cat %s | grep '!    total energy' >> energy-ecutwfc.data" % out_f_name)

ecut_wfc_all = [ ecut_wfc_min + i * ecut_wfc_step for i in range(n_test + 1)]
energy_all = []
with open("energy-ecutwfc.data", 'r') as fin:
    for line in fin:
        energy_all.append(float(line.split()[4]))

import matplotlib.pyplot as plt

plt.plot(ecut_wfc_all, energy_all)
plt.show()
