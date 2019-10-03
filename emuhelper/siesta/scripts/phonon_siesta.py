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
    python phonon_siesta.py xxx.xyz 
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
fc_displace = 0.04 # default value
fc_first = 1 # default
fc_last = fc_first # default

xyz = siesta_xyz(sys.argv[1])

system_name = "Phonon Calculation"
system_label = "phonon_calc"


if os.path.exists("./tmp-phonon"):
    shutil.rmtree("./tmp-phonon")
os.mkdir("./tmp-phonon")
os.chdir("./tmp-phonon")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
os.system("cp ../*.psf ./")

fdf_name = "phonon-calc.fdf"
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
# Phonon Calculation SETTING
with open(fdf_name, 'a') as fout:
    fout.write("MD.TypeOfRun FC\n")
    fout.write("MD.FCDispl %f Bohr\n" % fc_displace)
    fout.write("MD.FCFirst %d\n" % fc_first)
    fout.write("MD.FCLast %d\n" % fc_last)
    fout.write("\n")
    
# run the simulation
fdf_name = "phonon-calc.fdf"
out_f_name = "phonon-calc.out"
os.system("siesta < %s | tee %s" % (fdf_name, out_f_name))


# analyse the result
import matplotlib.pyplot as plt
