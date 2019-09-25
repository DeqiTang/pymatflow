#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import re
import pymatgen as mg

from emuhelper.qe.base.xyz import qe_xyz

"""
Usage:
    python phonon_with_phx_pwscf.py xxx.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential
    file for all the elements in the system is in the directory.

Note:
"""

# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])

ecut_wfc = 55
ecut_rho = ecut_wfc * 4

xyz = qe_xyz(sys.argv[1])

title = "Phonon with Ph.x"

base_prefix = "bfo"

if os.path.exists("./tmp-phonon-with-phx"):
    shutil.rmtree("./tmp-phonon-with-phx")
os.mkdir("./tmp-phonon-with-phx")
os.chdir("./tmp-phonon-with-phx")

os.system("cp ../*.UPF ./")

inp_name = "phonon.in"
with open(inp_name, 'w') as fout:
    fout.write("&control\n")
    fout.write("calculation = 'scf'\n")
    fout.write("title = '%s'\n" % title)
    fout.write("prefix = '%s'\n" % (base_prefix))
    fout.write("restart_mode = 'from_scratch'\n")
    fout.write("nstep = 300\n")
    fout.write("outdir = '%s'\n" % ("./tmp"))
    fout.write("pseudo_dir = './'\n")
    fout.write("wf_collect = .true.\n")
    fout.write("tstress = .true.\n")
    fout.write("tprnfor = .true.\n")
    fout.write("/\n")
    fout.write("\n")

    fout.write("&system\n")
    fout.write("ibrav = 0\n")
    fout.write("nat = %d\n" % xyz.natom)
    fout.write("ntyp = %d\n" % xyz.nspecies)
    fout.write("nspin = 1\n")
    #fout.write("nspin = 2\n")
    #fout.write("starting_magnetization(1) = 1\n")
    #fout.write("starting_magnetization(2) = 1\n")
    fout.write("ecutwfc = %d\n" % ecut_wfc)
    fout.write("ecutrho = %d\n" % ecut_rho)
    fout.write("input_DFT = 'PBE'\n")
    fout.write("occupations = 'smearing'\n")
    fout.write("degauss = 1.0d-4\n")
    fout.write("smearing = 'marzari-vanderbilt'\n")
    fout.write("/\n")
    fout.write("\n")

    fout.write("&electrons\n")
    fout.write("electron_maxstep = 300\n")
    #fout.write("conv_thr = 1.0d-10\n")
    fout.write("conv_thr = 1.0d-5\n")
    fout.write("mixing_mode = 'plain'\n")
    fout.write("mixing_beta = 0.3d0\n")
    fout.write("scf_must_converge = .true.\n")
    fout.write("/\n")
    fout.write("\n")

    fout.write("K_POINTS automatic\n")
    fout.write("1 1 1 0 0 0\n")
    fout.write("\n")
xyz.to_pwscf(inp_name)

