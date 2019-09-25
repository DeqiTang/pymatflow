#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.siesta.base.xyz import siesta_xyz



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


class siesta_xyz_phonopy(siesta_xyz):
    """
    """
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)

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


       
        

# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])
meshcutoff = 60
supercell_n = "1 1 2"
xyz = siesta_xyz_phonopy()

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

