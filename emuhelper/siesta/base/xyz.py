#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.base.atom import Atom
from emuhelper.base.xyz import base_xyz


"""
Usage:
    python converge_ecut.py xxx.xyz ecut_min ecut_max ecut_step
    xxx.xyz is the input structure file

    make sure the xyz structure file and the corresponding pseudopotential
    file for the elements of the system is in the directory.
"""


class siesta_xyz(base_xyz):
    """
    a representation of xyz file
    """
    def __init__(self, xyz_f):
        super().__init__(xyz_f)

    def to_fdf(self, fname):
        cell = self.cell
        with open(fname, 'a') as fout:
            fout.write("%block ChemicalSpeciesLabel\n")
            for element in self.specie_labels:
                fout.write("\t%d\t%d\t%s\n" % (self.specie_labels[element], mg.Element(element).number, element))
            fout.write("%endblock ChemicalSpeciesLabel\n")
            fout.write("\n")

            fout.write("%block PAO.BasisSizes\n")
            for element in self.specie_labels:
                fout.write("\t%s\tDZP\n" % element)
            fout.write("%endblock PAO.BasisSizes\n")
            fout.write("\n")

            fout.write("AtomicCoordinatesFormat ScaledCartesian\n")
            # 这里可以用ScaledCartesian也可以用Ang, 因为我的LatticeConstant 设置为1Ang
            # 这样ScaledCartesian以LatticeConstant扩展后的值实际上与Ang是一样的
            fout.write("AtomCoorFormatOut Ang\n")
            fout.write("LatticeConstant 1.00000 Ang\n")
            fout.write("\n")
            
            fout.write("%block LatticeVectors\n")
            fout.write("%f %f %f\n" % (cell[0], cell[1], cell[2]))
            fout.write("%f %f %f\n" % (cell[3], cell[4], cell[5]))
            fout.write("%f %f %f\n" % (cell[6], cell[7], cell[8]))
            fout.write("%endblock LatticeVectors\n")
            fout.write("\n")

            fout.write("%block AtomicCoordinatesAndAtomicSpecies\n")
            for atom in self.atoms:
                fout.write("%f\t%f\t%f\t" % (atom.x, atom.y, atom.z))
                fout.write(str(self.specie_labels[atom.name]))
                fout.write("\n")
            fout.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
            fout.write("\n")

       
