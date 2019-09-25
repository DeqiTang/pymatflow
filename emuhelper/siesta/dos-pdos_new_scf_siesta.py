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
    python dos_siesta.py xxx.xyz dos_emin dos_emax dos_peak_width dos_npoint
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
dos_emin = float(sys.argv[2]) 
dos_emax = float(sys.argv[3])
dos_peak_width = float(sys.argv[4])
dos_npoint = int(sys.argv[5])
xyz = siesta_xyz(sys.argv[1])

system_name = "Calculate DOS"
system_label = "dos"


if os.path.exists("./tmp-dos-pdos"):
    shutil.rmtree("./tmp-dos-pdos")
os.mkdir("./tmp-dos-pdos")
os.chdir("./tmp-dos-pdos")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../Na.psf", "Na.psf")
#shutil.copyfile("../K.psf", "K.psf")
#shutil.copyfile("../Nb.psf", "Nb.psf")
#shutil.copyfile("../O.psf", "O.psf")
os.system("cp ../*.psf ./")

fdf_name = "dos.fdf"
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
    fout.write("%block ProjectedDensityOfStates\n")
    fout.write("%f %f %f %d eV\n" % (dos_emin, dos_emax, dos_peak_width, dos_npoint))
    fout.write("%endblock ProjectedDensityOfStates\n")
    fout.write("%block PDOS.kgrid.MonkhorstPack\n")
    fout.write("10 0 0 0.5\n")
    fout.write("0 10 0 0.5\n")
    fout.write("0 0 10 0.5\n")
    fout.write("%endblock PDOS.kgrid.MonkhorstPack\n")
    fout.write("\n")

# run the simulation
fdf_name = "dos.fdf"
out_f_name = "dos.out"
os.system("siesta < %s > %s" % (fdf_name, out_f_name))


# analyse the result
import matplotlib.pyplot as plt
# plot dos
energy = []
states = []
with open(system_label+".DOS", 'r') as fin:
    for line in fin:
        energy.append(float(line.split()[0]))
        states.append(float(line.split()[1]))
plt.plot(energy, states)
plt.show()
