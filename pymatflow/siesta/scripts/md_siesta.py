#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from emuhelper.siesta.base.xyz import siesta_xyz


"""
Usage:
    python md_siesta.py xxx.xyz 
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

md_type = "Verlet" # Verlet, Nose, ParrinelloRahman, NoseParrinelloRahman, Anneal
max_step = 20 # for MD.Steps: maximum step
md_final_t_step = max_step # default is MD.Steps 
xyz = siesta_xyz(sys.argv[1])

system_name = "Molecular Dynamics"
system_label = "md"


if os.path.exists("./tmp-md"):
    shutil.rmtree("./tmp-md")
os.mkdir("./tmp-md")
os.chdir("./tmp-md")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
os.system("cp ../*.psf ./")


fdf_name = "md.fdf"
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
# MD SETTING
with open(fdf_name, 'a') as fout:
    fout.write("MD.TypeOfRun %s\n" % md_type)
    fout.write("MD.InitialTimeStep 1\n")
    fout.write("MD.FinalTimeStep %d\n" % md_final_t_step)
    fout.write("MD.LengthTimeStep 1 fs\n")
    fout.write("MD.InitialTemperature 0 K\n")
    fout.write("MD.TargetTemperature 0 K\n")
    
# run the simulation
fdf_name = "md.fdf"
out_f_name = "md.out"
os.system("siesta < %s | tee %s" % (fdf_name, out_f_name))


# analyse the result
import matplotlib.pyplot as plt
