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
    python chem_analysis_siesta.py xxx.xyz 
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

system_name = "Chemistry Analysis"
system_label = "chem"


if os.path.exists("./tmp-chem-analysis"):
    shutil.rmtree("./tmp-chem-analysis")
os.mkdir("./tmp-chem-analysis")
os.chdir("./tmp-chem-analysis")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
os.system("cp ../*.psf ./")

fdf_name = "chem-analysis.fdf"
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

# run the simulation
fdf_name = "chem-analysis.fdf"
out_f_name = "chem-analysis.out"
os.system("siesta < %s > %s" % (fdf_name, out_f_name))


# analyse the result
import matplotlib.pyplot as plt

# Analysis the result of: Mulliken charges and overlap populations

# Analysis the result of: Voronoi and Hirshfeld atomic population analysis

# Analysis the result of: Crystal-Orbital overlap and hamilton populations(COOP/COHP)


