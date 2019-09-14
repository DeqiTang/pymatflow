#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
    python go_siesta.py xxx.xyz 
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
            fout.write("4.000175 0 0\n")
            fout.write("0 4.000175 0\n")
            fout.write("0 0 4.000175\n")
            fout.write("%endblock LatticeVectors\n")
            fout.write("\n")

            fout.write("%block AtomicCoordinatesAndAtomicSpecies\n")
            for atom in self.atoms:
                fout.write("%f\t%f\t%f\t" % (atom.x, atom.y, atom.z))
                fout.write(str(xyz.specie_labels[atom.name]))
                fout.write("\n")
            fout.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
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
meshcutoff = 60

geo_opt_type = "CG" # CG, Broyden, 
if_vc = "true"
if_cv = "false"
max_force = 0.01 # default value 0.04
max_stress = 0.1 # default value 1
max_step = 20
xyz = XYZ()

system_name = "Geometric Optimization"
system_label = "geo_opt"


if os.path.exists("./tmp-geo-opt"):
    shutil.rmtree("./tmp-geo-opt")
os.mkdir("./tmp-geo-opt")
os.chdir("./tmp-geo-opt")
shutil.copyfile("../H.psf", "H.psf")
shutil.copyfile("../Li.psf", "Li.psf")

fdf_name = "geo-opt.fdf"
with open(fdf_name, 'w') as fout:
    fout.write("SystemName\t%s\n" % (system_name))
    fout.write("SystemLabel\t%s\n" % (system_label))
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
# GEO OPT SETTING
with open(fdf_name, 'a') as fout:
    fout.write("MD.TypeOfRun %s\n" % geo_opt_type)
    fout.write("MD.VariableCell %s\n" % if_vc)
    fout.write("MD.ConstantVolume %s\n" % if_cv)
    fout.write("MD.MaxForceTol %f eV/Ang\n" % max_force)
    fout.write("MD.MaxStressTol %f GPa\n" % max_stress)
    fout.write("MD.Steps %d\n" % max_step)
    fout.write("MD.MaxDispl 0.2 Bohr # default value\n")
    fout.write("MD.PreconditionVariableCell 5 Ang # default value\n")
    
# run the simulation
fdf_name = "geo-opt.fdf"
out_f_name = "geo-opt.out"
os.system("siesta < %s > %s" % (fdf_name, out_f_name))


# analyse the result
import matplotlib.pyplot as plt