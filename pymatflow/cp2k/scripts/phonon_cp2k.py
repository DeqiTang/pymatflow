#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.cp2k.base.xyz import cp2k_xyz


"""
Usage:
    python phonon_cp2k.py xxx.xyz 
    xxx.xyz is the input structure file

    make sure the xyz structure file and pseudopotential file
    for all the elements of the system is in the directory.

Dependencies:
    pip install --user phonopy
    pip install --user cp2k_tools

Note:
    phonopy read the xxx.inp and it can only read the system structure
    by COORD specified in SUBSYS. So I can not use TOPOLOGY.
    PLUS: only scaled coordinates are currently supported!

References:
    https://www.cp2k.org/exercises:2018_uzh_cmest:phonon_calculation
"""


class cp2k_xyz_phonopy(cp2k_xyz):
    """
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def to_subsys_phonopy(self, fname):
        cell = self.cell
        with open(fname, 'a') as fout:
            fout.write("\t&SUBSYS\n")
            for element in self.specie_labels:
                fout.write("\t\t&KIND %s\n" % element)
                fout.write("\t\t\tBASIS_SET DZVP-MOLOPT-SR-GTH\n")
                fout.write("\t\t\tPOTENTIAL GTH-PBE\n")
                fout.write("\t\t&END KIND\n")
            fout.write("\t\t&CELL\n")
            fout.write("\t\t\tABC %f %f %f\n" % (cell[0], cell[4], cell[8]))
            #fout.write("\t\t\tA %f %f %f\n" % (cell[0], cell[1], cell[2]))
            #fout.write("\t\t\tB %f %f %f\n" % (cell[3], cell[4], cell[5]))
            #fout.write("\t\t\tC %f %f %f\n" % (cell[6], cell[7], cell[8]))
            fout.write("\t\t&END CELL\n")
            #fout.write("\t\t&TOPOLOGY\n")
            #fout.write("\t\t\tCOORD_FILE_FORMAT xyz\n")
            #fout.write("\t\t\tCOORD_FILE_NAME %s\n" % sys.argv[1])
            #fout.write("\t\t&END TOPOLOGY\n")
            fout.write("\t\t&COORD\n")
            fout.write("\t\t\tSCALED .TRUE.\n")
            for atom in self.atoms:
                fout.write("\t\t\t%s\t%f\t%f\t%f\n" % (atom.name, atom.x/cell[0], atom.y/cell[4], atom.z/cell[8]))
            fout.write("\t\t&END COORD\n")
            fout.write("\t&END SUBSYS\n")
            fout.write("\n")

    def print_kinds(self, fname):
        with open(fname, 'a') as fout:
            for element in self.specie_labels:
                fout.write("\t\t&KIND %s\n" % element)
                fout.write("\t\t\tBASIS_SET DZVP-MOLOPT-SR-GTH\n")
                fout.write("\t\t\tPOTENTIAL GTH-PBE\n")
                fout.write("\t\t&END KIND\n")


# OK now we can use XYZ class to extract information 
# from the xyz file: sys.argv[1]

#cutoff_min = int(sys.argv[2]) # in Ry: 1 Ry = 13.6 ev
#cutoff_max = int(sys.argv[3])
#cutoff_step = int(sys.argv[4])
#rel_cutoff = int(sys.argv[5])
cutoff = 60
rel_cutoff = 30

supercell_n = "1 1 1"

xyz = cp2k_xyz_phonopy(sys.argv[1])

base_project_name = "phonon"

if os.path.exists("./tmp-phonon"):
    shutil.rmtree("./tmp-phonon")
os.mkdir("./tmp-phonon")
os.chdir("./tmp-phonon")
shutil.copyfile("../%s" % sys.argv[1], "%s" % sys.argv[1])
#shutil.copyfile("../Li.psf", "Li.psf")

inp_name = "phonon.inp"
with open(inp_name, 'w') as fout:
    fout.write("&GLOBAL\n")
    fout.write("\tPROJECT\t%s\n" % (base_project_name))
    fout.write("\tRUN_TYPE ENERGY_FORCE\n")
    fout.write("\tPRINT_LEVEL LOW\n")
    fout.write("&END GLOBAL\n")
    fout.write("\n")

    fout.write("&FORCE_EVAL\n")
    fout.write("\tMETHOD Quickstep\n")
# subsys
xyz.to_subsys_phonopy(inp_name)
# end subsys
with open(inp_name, 'a') as fout:
    # dft
    fout.write("\t&dft\n")
    fout.write("\t\tBASIS_SET_FILE_NAME BASIS_MOLOPT\n")
    fout.write("\t\tPOTENTIAL_FILE_NAME GTH_POTENTIALS\n")
    fout.write("\t\t&QS\n")
    fout.write("\t\t\tEPS_DEFAULT 1.0E-10\n")
    fout.write("\t\t&END QS\n")
    fout.write("\t\t&MGRID\n")
    fout.write("\t\t\tNGRIDS 4\n")
    fout.write("\t\t\tCUTOFF %d\n" % cutoff)
    fout.write("\t\t\tREL_CUTOFF %d\n" % rel_cutoff)
    fout.write("\t\t&END MGRID\n")
    fout.write("\t\t&XC\n")
    fout.write("\t\t\t&XC_FUNCTIONAL PBE\n")
    fout.write("\t\t\t&END XC_FUNCTIONAL\n")
    fout.write("\t\t&END XC\n")
    fout.write("\t\t&SCF\n")
    fout.write("\t\t\tSCF_GUESS ATOMIC\n")
    fout.write("\t\t\tEPS_SCF 1.0E-06\n")
    fout.write("\t\t\tMAX_SCF 200\n")
    fout.write("\t\t\tADDED_MOS 100\n")
    fout.write("\t\t\t&DIAGONALIZATION ON\n")
    fout.write("\t\t\t\tALGORITHM STANDARD\n")
    fout.write("\t\t\t&END DIAGONALIZATION\n")
    fout.write("\t\t\t&SMEAR ON\n")
    fout.write("\t\t\t\tMETHOD FERMI_DIRAC\n")
    fout.write("\t\t\t\tELECTRONIC_TEMPERATURE [K] 300\n")
    fout.write("\t\t\t&END SMEAR\n")
    fout.write("\t\t\t&MIXING T\n")
    fout.write("\t\t\t\tMETHOD BROYDEN_MIXING\n")
    fout.write("\t\t\t\tALPHA 0.4\n")
    fout.write("\t\t\t\tNBROYDEN 8\n")
    fout.write("\t\t\t&END MIXING\n")
    fout.write("\t\t&END SCF\n")
    fout.write("\t&END DFT\n")
    # end dft
    fout.write("&END FORCE_EVAL\n")


# run the simulation
os.system("phonopy --cp2k -c %s -d --dim='%s'" % (inp_name, supercell_n))
# now supercell-00x.inp is generated which will be used to construct input for cp2k
os.system("ls | grep 'supercell-' > geo.data")
disps = []
with open("geo.data", 'r') as fin:
    for line in fin:
        disps.append(line.split(".")[0].split("-")[1])

for disp in disps:
    in_name = "supercell-%s.inp" % disp
    if os.path.exists(in_name) is not True:
        break
    tmp_file = "supercell-%s.tmp.txt" % disp
    shutil.copyfile(in_name, tmp_file)
    with open(in_name, 'w') as fout:
        fout.write("&GLOBAL\n")
        fout.write("\tPROJECT\t%s\n" % (base_project_name+"-supercell-"+disp))
        fout.write("\tRUN_TYPE ENERGY_FORCE\n")
        fout.write("\tPRINT_LEVEL LOW\n")
        fout.write("&END GLOBAL\n")
        fout.write("\n")

        fout.write("&FORCE_EVAL\n")
        fout.write("\tMETHOD Quickstep\n")
        fout.write("\t&SUBSYS\n")
    xyz.print_kinds(in_name)
    os.system("cat %s | sed '1d;2d;3d;4d;5d;6d;7d' | sed '$d' | sed '$d' | sed '$d' | sed '$d' | sed '$d' >> %s" % (tmp_file, in_name))
    with open(in_name, 'a') as fout:
        fout.write("\t&END SUBSYS\n")
        # dft
        fout.write("\t&DFT\n")
        fout.write("\t\tBASIS_SET_FILE_NAME BASIS_MOLOPT\n")
        fout.write("\t\tPOTENTIAL_FILE_NAME GTH_POTENTIALS\n")
        fout.write("\t\t&QS\n")
        fout.write("\t\t\tEPS_DEFAULT 1.0E-10\n")
        fout.write("\t\t&END QS\n")
        fout.write("\t\t&MGRID\n")
        fout.write("\t\t\tNGRIDS 4\n")
        fout.write("\t\t\tCUTOFF %d\n" % cutoff)
        fout.write("\t\t\tREL_CUTOFF %d\n" % rel_cutoff)
        fout.write("\t\t&END MGRID\n")
        fout.write("\t\t&XC\n")
        fout.write("\t\t\t&XC_FUNCTIONAL PBE\n")
        fout.write("\t\t\t&END XC_FUNCTIONAL\n")
        fout.write("\t\t&END XC\n")
        fout.write("\t\t&SCF\n")
        fout.write("\t\t\tSCF_GUESS ATOMIC\n")
        fout.write("\t\t\tEPS_SCF 1.0E-06\n")
        fout.write("\t\t\tMAX_SCF 200\n")
        fout.write("\t\t\tADDED_MOS 100\n")
        fout.write("\t\t\t&DIAGONALIZATION ON\n")
        fout.write("\t\t\t\tALGORITHM STANDARD\n")
        fout.write("\t\t\t&END DIAGONALIZATION\n")
        fout.write("\t\t\t&SMEAR ON\n")
        fout.write("\t\t\t\tMETHOD FERMI_DIRAC\n")
        fout.write("\t\t\t\tELECTRONIC_TEMPERATURE [K] 300\n")
        fout.write("\t\t\t&END SMEAR\n")
        fout.write("\t\t\t&MIXING T\n")
        fout.write("\t\t\t\tMETHOD BROYDEN_MIXING\n")
        fout.write("\t\t\t\tALPHA 0.4\n")
        fout.write("\t\t\t\tNBROYDEN 8\n")
        fout.write("\t\t\t&END MIXING\n")
        fout.write("\t\t&END SCF\n")
        fout.write("\t&END DFT\n")
        # end dft
        fout.write("\t&PRINT\n")
        fout.write("\t\t&FORCES\n")
        fout.write("\t\t\tFILENAME forces\n")
        fout.write("\t\t&END FORCES\n")
        fout.write("\t&END PRINT\n")
        fout.write("&END FORCE_EVAL\n")

for disp in disps:
    in_name = "supercell-%s.inp" % disp
    if os.path.exists(in_name) is not True:
        break
    os.system("cp2k.psmp -in %s | tee %s" % (in_name, in_name+".out"))


phonopy_command = "phonopy --cp2k -f "
for disp in disps:
    f_name = base_project_name + "-supercell-" + disp + "-forces-1_0.xyz"
    if os.path.exists(f_name) is not True:
        break
    phonopy_command = phonopy_command + f_name + " "

os.system(phonopy_command)

# get the band structure
# 注意--pa设置Primitive Axis要设置正确! --band 控制了声子谱的图示
os.system("phonopy --cp2k -c %s -p --dim='%s' --pa='1 0 0 0 1 0 0 0 1' --band='1/2 1/2 1/2 0 0 0 1/2 0 1/2'" % (inp_name, supercell_n))

# analyse the result

import matplotlib.pyplot as plt

