#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.abinit.base.xyz import abinit_xyz

"""
Usage:
    python neb_abinit.py initial.xyz final.xyz
    xxx.xyz is the input structure file

    make sure the xyz structure file and the pseudopotential file
    for all the elements of the system is in the directory.
Note:
    参考: https://docs.abinit.org/topics/TransPath/
"""

class abinit_xyz_neb(abinit_xyz):
    """

    """
    def __init__(self, xyz_f):
        super().__init__(xyz_f)
    #def __init__(self, *args, **kwargs):
    #    super().__init__(*args, **kwargs)

    def to_xangst_lastimg(self, fname):
        with open(fname, 'a') as fout:
            fout.write("xangst_lastimg\n")
            for atom in self.atoms:
                fout.write("%f %f %f\n" % (atom.x, atom.y, atom.z))
            fout.write("\n")


# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#ecut_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#ecut_max = int(sys.argv[3])
#ecut_step = int(sys.argv[4])
cutoff = 40

xyz_initial = abinit_xyz_neb(sys.argv[1])
xyz_final = abinit_xyz_neb(sys.argv[2])

base_project_name = "neb-transpath"
if os.path.exists("./tmp-neb"):
    shutil.rmtree("./tmp-neb")
os.mkdir("./tmp-neb")
os.chdir("./tmp-neb")
#shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")
#shutil.copyfile("../Li.psp8", "Li.psp8")
#shutil.copyfile("../H.psp8", "H.psp8")
os.system("cp ../*.psp8 ./")


inp_name = "neb.in"
files_name = "neb.files"
with open(files_name, 'w') as fout:
    fout.write(inp_name)
    fout.write("\n")
    fout.write("neb.out\n")
    fout.write("nebi\n")
    fout.write("nebo\n")
    fout.write("temp\n")
    for element in xyz_initial.specie_labels:
        fout.write("%s\n" % (element + ".psp8"))
    #
with open(inp_name, 'w') as fout:
    fout.write("ecut %d\n" % cutoff)
    fout.write("kptopt 1\n")
    fout.write("ngkpt 1 1 1\n")
    fout.write("occopt 3\n") # fermi dirac smearing of occupation
    fout.write("nstep 100\n")
    fout.write("toldff 1.0d-6\n")
    fout.write("nband 10\n")
    fout.write("diemac 2.0\n")
    fout.write("nimage 12\n")
    fout.write("imgmov 5\n")
    fout.write("ntimimage 50\n")
    fout.write("tolimg 0.0001\n") # Tol. criterion (will stop when average energy of cells < tolimg)
    fout.write("dynimage 0 10*1 0\n")  # Keep first and last images fixed
    fout.write("fxcartfactor 1.0\n")
    fout.write("prtvolimg 0\n")
    fout.write("\n")
xyz_initial.to_abinit(inp_name)
xyz_final.to_xangst_lastimg(inp_name)


# run the simulation
#out_f_name = "geo-opt-calc.out.log"
os.system("abinit < %s" % (files_name))

# analyse the result

import matplotlib.pyplot as plt

