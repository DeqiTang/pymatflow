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
    python charge-density_potentials_siesta.py xxx.xyz 
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

system_name = "Output of charge densities and potentials on the grid"
system_label = "charge-density-potentials"


if os.path.exists("./tmp-charge-density-potentials"):
    shutil.rmtree("./tmp-charge-density-potentials")
os.mkdir("./tmp-charge-density-potentials")
os.chdir("./tmp-charge-density-potentials")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
os.system("cp ../*.psf ./")

fdf_name = "charge-density-potentials.fdf"
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
    # Output of charge densities and potentials on the grid
    fout.write("SaveRho true\n")
    fout.write("SaveDeltaRho true\n")
    fout.write("SaveRhoXC true\n")
    fout.write("SaveElectricalstaticPotential true\n")
    fout.write("SaveNeutralAtomPotential true\n")
    fout.write("SaveTotalPotential true\n")
    fout.write("SaveIonicCharge true\n")
    fout.write("SaveTotalCharge true\n")
    fout.write("SaveBaderCharge true\n")
    fout.write("AnalyzeChargeDensityOnly false\n")
    fout.write("SaveInitialChargeDensity false\n")
    fout.write("\n")

# run the simulation
fdf_name = "charge-density-potentials.fdf"
out_f_name = "charge-density-potentials.out"
os.system("siesta < %s | tee %s" % (fdf_name, out_f_name))


# analyse the result
import matplotlib.pyplot as plt

# Analysis the result: Output of charge densities and potentials on the grid

# Analysis the result: Bader charge

# Analysis the result: 
