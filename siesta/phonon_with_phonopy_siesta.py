#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
Usage:
    python phonon_with_phonopy_siesta.py xxx.xyz 
    xxx.xyz is the input structure file

    make sure the xyz structure file and the corresponding pseudopotential
    file for the elements of the system is in the directory.
Note:
    参考:
    https://atztogo.github.io/phonopy/siesta.html
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
            # 这里可以用ScaledCartesian也可以用Ang, 因为我的LatticeConstant 设置为1Ang
            # 这样ScaledCartesian以LatticeConstant扩展后的值实际上与Ang是一样的
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
    def for_phonopy(self, fname):
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
meshcutoff = 60
supercell_n = "1 1 2"
xyz = XYZ()

system_name = "Phonon spectroscopy with phonopy"
system_label = "phonon"


if os.path.exists("./tmp-phonon-with-phonopy"):
    shutil.rmtree("./tmp-phonon-with-phonopy")
os.mkdir("./tmp-phonon-with-phonopy")
os.chdir("./tmp-phonon-with-phonopy")
#shutil.copyfile("../Fe.psf", "Fe.psf")
#shutil.copyfile("../Bi.psf", "Bi.psf")
#shutil.copyfile("../O.psf", "O.psf")
os.system("cp ../*.psf ./")

head_fdf_name = "head.fdf"
with open(head_fdf_name, 'w') as fout:
    fout.write("SystemName\t%s\n" % (system_name))
    fout.write("SystemLabel\t%s\n" % (system_label))
    fout.write("NumberOfSpecies\t%d\n" %  xyz.nspecies)
    #fout.write("NumberOfAtoms\t%d\n" % xyz.natom)        
    fout.write("\n")
xyz.for_phonopy(head_fdf_name)
with open(head_fdf_name, 'a') as fout:
    #---------------------------------
    # DFT settings:
    # --------------------------------
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
    # Mulliken charges and overlap populations
    fout.write("WriteMullikenPop 3\n")
    fout.write("MullikenInSCF false\n")
    fout.write("SpinInSCF true\n")
    # Voronoi and Hirshfeld atomic population analysis
    fout.write("WriteHirshfeldPop true\n")
    fout.write("WriteVoronoiPop true\n")
    fout.write("PartialChargesAtEveryGeometry false\n")
    fout.write("PartialChargesAtEverySCFStep false\n")
    # Crystal-Orbital overlap and hamilton populations(COOP/COHP)
    fout.write("COOP.Write true\n")
    fout.write("#WFS.Energy.Min")   # default −∞
    fout.write("#WFS.Energy.Max")   # default ∞
    fout.write("\n")

pos_fdf_name = "pos.fdf"
xyz.to_fdf(pos_fdf_name)

# set up the Phonopy calculation
os.system("phonopy --siesta -d --dim='%s' -c %s" % (supercell_n, pos_fdf_name))
os.system("ls | grep 'supercell-' > pos.data")
disp_dirs = []
with open("pos.data", 'r') as fin:
    for line in fin:
        disp_dirs.append(line.split(".")[0].split("-")[1])
for disp in disp_dirs:
    os.mkdir("disp-%s" % disp)
    os.system("cat %s supercell-%s.fdf > ./disp-%s/supercell-%s.fdf" % (head_fdf_name, disp, disp, disp))
    os.system("cp *.psf ./disp-%s/" % disp)
    os.system("rm supercell-%s.fdf" % disp)
# run every disp
for disp in disp_dirs:
    os.chdir("disp-%s" % disp)
    os.system("siesta < supercell-%s.fdf > supercell-%s.out" % (disp, disp))
    os.chdir("../")

# analyse the result
import matplotlib.pyplot as plt
# create FORCE_SETS
os.system("phonopy --siesta -f disp-{001..%s}/%s.FA" % (disp_dirs[-1], system_label))
# plot the phonon band
with open("band.conf", 'w') as fout:
    fout.write("ATOM_NAME =")
    for element in xyz.specie_labels:
        fout.write(" %s" % element)
    fout.write("\n")
    fout.write("DIM = %s\n" % (supercell_n))
    fout.write("BAND = 0.5 0.5 0.5 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.5 0.0\n")
os.system("phonopy --siesta -c %s -p band.conf" % pos_fdf_name)

