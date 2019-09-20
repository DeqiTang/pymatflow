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

    make sure the xyz structure file and the corresponding pseudopotential
    file for the elements of the system is in the directory.
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

    def to_fdf(self, fname):
        cell = self.cell
        with open(fname, 'a') as fout:
            fout.write("%block ChemicalSpeciesLabel\n")
            for element in self.specie_labels:
                fout.write("\t%d\t%d\t%s\n" % (self.specie_labels[element], mg.Element(element).number, element))
            fout.write("%endblock ChemicalSpeciesLabel\n")
            fout.write("\n")

            fout.write("%block PAO.BasisSizes\n")
            for element in self.specie_labels:
                fout.write("\t%s\tDZP\n" % element)
            fout.write("%endblock PAO.BasisSizes\n")
            fout.write("\n")

            fout.write("AtomicCoordinatesFormat ScaledCartesian\n")
            fout.write("AtomCoorFormatOut Ang\n")
            fout.write("LatticeConstant 1.00000 Ang\n")
            fout.write("\n")
            
            fout.write("%block LatticeVectors\n")
            fout.write("%f %f %f\n" % (cell[0], cell[1], cell[2]))
            fout.write("%f %f %f\n" % (cell[3], cell[4], cell[5]))
            fout.write("%f %f %f\n" % (cell[6], cell[7], cell[8]))
            fout.write("%endblock LatticeVectors\n")
            fout.write("\n")

            fout.write("%block AtomicCoordinatesAndAtomicSpecies\n")
            for atom in self.atoms:
                fout.write("%f\t%f\t%f\t" % (atom.x, atom.y, atom.z))
                fout.write(str(xyz.specie_labels[atom.name]))
                fout.write("\n")
            fout.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
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

ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
ecut_max = int(sys.argv[3])
ecut_step = int(sys.argv[4])


xyz = XYZ()

system_name = "Test Meshcutoff Value"
system_label = "TestEcut"


if os.path.exists("./tmp-ecut"):
    shutil.rmtree("./tmp-ecut")
os.mkdir("./tmp-ecut")
os.chdir("./tmp-ecut")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../C.psf", "C.psf")
shutil.copyfile("../Na.psf", "Na.psf")
shutil.copyfile("../K.psf", "K.psf")
shutil.copyfile("../Nb.psf", "Nb.psf")
shutil.copyfile("../O.psf", "O.psf")


n_test = int((ecut_max - ecut_min) / ecut_step)
for i in range(n_test + 1):
    meshcutoff = int(ecut_min + i * ecut_step)
    fdf_name = "test-ecut-%d.fdf" % meshcutoff
    with open(fdf_name, 'w') as fout:
        fout.write("SystemName\t%s\n" % (system_name + str(meshcutoff)))
        fout.write("SystemLabel\t%s\n" % (system_label + str(meshcutoff)))
        fout.write("NumberOfSpecies\t%d\n" %  xyz.nspecies)
        fout.write("NumberOfAtoms\t%d\n" % xyz.natom)        
        fout.write("\n")
    xyz.to_fdf(fdf_name)
    #---------------------------------
    # DFT settings:
    # --------------------------------
    with open(fdf_name, 'a') as fout:
        fout.write("# DFT settings\n")
        fout.write("\n")
        fout.write("XC.functional  GGA\n")
        fout.write("XC.authors  PBE\n")
        fout.write("DM.Tolerance  1.d-6\n")
        fout.write("MaxSCFIterations  300\n")
        fout.write("DM.MixingWeight  0.1\n")
        fout.write("DM.NumberPulay  5\n")
        fout.write("DM.AllowExtrapolation  true\n")
        # do not DM.UseSaveDM so that different ecut does not affect each other
        fout.write("DM.UseSaveDM  false\n")         
        fout.write("SolutionMethod diagon\n")
        fout.write("MeshCutoff  %d Ry\n" % meshcutoff)
        fout.write("\n")

# run the simulation
for i in range(n_test + 1):
    meshcutoff = int(ecut_min + i * ecut_step)
    fdf_name = "test-ecut-%d.fdf" % meshcutoff
    out_f_name = "test-ecut-%d.out" % meshcutoff
    os.system("siesta < %s > %s" % (fdf_name, out_f_name))


# analyse the result
for i in range(n_test + 1):
    meshcutoff = int(ecut_min + i * ecut_step)
    out_f_name = "test-ecut-%d.out" % meshcutoff
    os.system("cat %s | grep 'Total =' >> energy-ecut.data" % out_f_name)

ecut = [ ecut_min + i * ecut_step for i in range(n_test + 1)]
energy = []
with open("energy-ecut.data", 'r') as fin:
    for line in fin:
        energy.append(float(line.split()[3]))

import matplotlib.pyplot as plt

plt.plot(ecut, energy)
plt.show()
