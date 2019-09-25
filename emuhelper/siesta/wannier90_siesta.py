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
    python wannier90_siesta.py xxx.xyz 
    xxx.xyz is the input structure file

    make sure the xyz structure file and the corresponding pseudopotential
    file for the elements of the system is in the directory.
"""



# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])
meshcutoff = 60
#dos_emin = float(sys.argv[2]) 
#dos_emax = float(sys.argv[3])
#dos_peak_width = float(sys.argv[4])
#dos_npoint = int(sys.argv[5])
xyz = siesta_xyz(sys.argv[1])

system_name = "Maximally Localized Wannier Functions"
system_label = "wannier90-functions"


if os.path.exists("./tmp-wannier90"):
    shutil.rmtree("./tmp-wannier90")
os.mkdir("./tmp-wannier90")
os.chdir("./tmp-wannier90")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
os.system("cp ../*.psf ./")

fdf_name = "wannier90-output.fdf"
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
    # wannier90: Maximally Localized Wannier Functions
    fout.write("Siesta2Wannier90.WriteMmm true\n")
    fout.write("Siesta2Wannier90.WriteAmm true\n")
    fout.write("Siesta2Wannier90.WriteEig true\n")
    fout.write("Siesta2Wannier90.WriteUnk true\n")
    fout.write("Siesta2Wannier90.UnkGrid1 10\n")
    fout.write("Siesta2Wannier90.UnkGrid2 10\n")
    fout.write("Siesta2Wannier90.UnkGrid3 10\n")
    fout.write("Siesta2Wannier90.UnkGridBinary true\n")
    fout.write("#Siesta2Wannier90.NumberOfBands # default: occupied bands\n") # default: occupied bands
    fout.write("#Siesta2Wannier90.NumberOfBandsUp # default: Siesta2Wannier90.NumberOfBands")
    fout.write("#Siesta2Wannier90.NumberOfBandsDown # default: Siesta2Wannier90.NumberOfBands")
    fout.write("\n")

# run the simulation
fdf_name = "wannier90-output.fdf"
out_f_name = "wannier90-output.out"
os.system("siesta < %s > %s" % (fdf_name, out_f_name))


# analyse the result
import matplotlib.pyplot as plt

