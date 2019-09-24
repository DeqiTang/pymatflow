#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from base.xyz import siesta_xyz



"""
Usage:
    python geo_opt_siesta.py xxx.xyz 
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

geo_opt_type = "CG" # CG, Broyden, 
if_vc = "true"
if_cv = "false"
max_force = 0.01 # default value 0.04
max_stress = 0.1 # default value 1
max_step = 20
xyz = siesta_xyz(sys.argv[1])

system_name = "Geometric Optimization"
system_label = "geo_opt"


if os.path.exists("./tmp-geo-opt"):
    shutil.rmtree("./tmp-geo-opt")
os.mkdir("./tmp-geo-opt")
os.chdir("./tmp-geo-opt")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../C.psf", "C.psf")
os.system("cp ../*.psf ./")

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
    # write the final structure to an xyz file
    fout.write("WriteCoorXmol true\n")
    fout.write("\n")
    
# run the simulation
fdf_name = "geo-opt.fdf"
out_f_name = "geo-opt.out"
os.system("siesta < %s > %s" % (fdf_name, out_f_name))



# analyse the result

import matplotlib.pyplot as plt

os.system("cat %s | grep 'siesta: E_KS(eV) =' > energy-per-geo-step.data" % (out_f_name))

energies = []
with open("energy-per-geo-step.data", 'r') as fin:
    for line in fin:
        energies.append(float(line.split()[3]))

steps = [i for i in range(len(energies))]
plt.plot(steps, energies)
plt.show()
