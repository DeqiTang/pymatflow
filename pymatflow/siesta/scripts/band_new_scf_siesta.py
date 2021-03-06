#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

from emuhelper.siesta.base.xyz import siesta_xyz

"""
Usage:
    python dos_siesta.py xxx.xyz dos_emin dos_emax dos_peak_width dos_npoint
    xxx.xyz is the input structure file

    make sure the xyz structure file and the corresponding pseudopotential
    file for the elements of the system is in the directory.
Notes:
    when you calculate the band you must physically understand your system.
    if your system defined is really a periodic system(with proper cell) you get solid bands like usually
    but if your defined system is not physically periodic(ie. use very large cell so that the system is
    actually isolated), you will not get bands but actually some energy levels!
    so if you are calculating the bands make sure you set the appropriate cell parameters.
    PS: SIESTA will use the cell you defined to duplicate the system periodically, but if your cell
    is significantly larger than your system, it is actually a isolated cluster or something physically
    so there is no bands but energy levels.
    And so as the case of DOS calculation
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

system_name = "Calculate bands"
system_label = "bands"


if os.path.exists("./tmp-bands"):
    shutil.rmtree("./tmp-bands")
os.mkdir("./tmp-bands")
os.chdir("./tmp-bands")
#shutil.copyfile("../H.psf", "H.psf")
#shutil.copyfile("../Li.psf", "Li.psf")
os.system("cp ../*.psf ./")

fdf_name = "bands.fdf"
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
    fout.write("BandLinesScale pi/a\n")
    fout.write("%block BandLines\n")
    fout.write("1 1.0 1.0 1.0 L # begin at L\n")
    fout.write("20 0.0 0.0 0.0 \Gamma # 20 points from L to gamma\n")
    fout.write("25 2.0 0.0 0.0 X # 25 points from gamma to X \n")
    fout.write("30 2.0 2.0 2.0 \Gamma # 30 points from X to gamma\n")
    fout.write("%endblock BandLines\n")
    fout.write("WriteKbands false\n") # do not write kbands to the standard out
    fout.write("WriteBands false\n") # do not write bands to the standard out
    fout.write("\n")

# run the simulation
fdf_name = "bands.fdf"
out_f_name = "bands.out"
os.system("siesta < %s | tee %s" % (fdf_name, out_f_name))


# analyse the result
import matplotlib.pyplot as plt

# plot band
os.system("gnubands < bands.bands > bands.data")
nk = 76
nbands = 40
kpoints = []
band_ens = []

with open("bands.data", 'r') as fin:
    for i in range(16):
        fin.readline()
    for i in range(nbands):
        band_ens.append([])
        if i == 0:
            for j in range(nk):
                line = fin.readline()
                kpoints.append(float(line.split()[0]))
                band_ens[i].append(float(line.split()[1]))
        else:
            for j in range(nk):
                line = fin.readline()
                band_ens[i].append(float(line.split()[1]))
        fin.readline()
        fin.readline()

for i in range(10):
    plt.plot(kpoints, band_ens[i][:])
plt.show()
