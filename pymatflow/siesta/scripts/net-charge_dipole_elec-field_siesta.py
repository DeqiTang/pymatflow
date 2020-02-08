#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from emuhelper.siesta.base.xyz import siesta_xyz



"""
Usage:
    python net-charge-dipole-elec-field.py xxx.xyz 
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

system_name = "Systems with net charge or dipole, and electric fields"
system_label = "net-charge-dipole-elec-field"


if os.path.exists("./tmp-net-charge-dipole-elec-field"):
    shutil.rmtree("./tmp-net-charge-dipole-elec-field")
os.mkdir("./tmp-net-charge-dipole-elec-field")
os.chdir("./tmp-net-charge-dipole-elec-field")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
os.system("cp ../*.psf ./")

fdf_name = "net-charge-dipole-elec-field.fdf"
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
    # Systems with net charge or dipole, and electric fields
    fout.write("NetCharge 1.0\n")
    fout.write("SimulateDoping true\n")
    fout.write("%block ExternalElectricField\n")
    fout.write("  0.000 0.000 0.50 V/Ang\n")
    fout.write("%endblock ExternalELectricField\n")
    fout.write("SlabDipoleCorrection true\n")
    #fout.write("%block Geometry.Hartree\n")
    #fout.write("%endblock Geometry.Hartree\n")
    #fout.write("%block Geometry.Charge\n")
    #fout.write("%endblock Geometry.Charge\n")
    fout.write("\n")

# run the simulation
fdf_name = "net-charge-dipole-elec-field.fdf"
out_f_name = "net-charge-dipole-elec-field.out"
os.system("siesta < %s | tee %s" % (fdf_name, out_f_name))


# analyse the result
import matplotlib.pyplot as plt

