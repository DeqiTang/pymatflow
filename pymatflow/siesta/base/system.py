#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from pymatflow.base.atom import Atom
from pymatflow.base.xyz import base_xyz


"""
Usage:
"""


class siesta_system:
    """
    """
    def __init__(self, xyz_f, name="Siesta Job", label="siesta"):
        self.xyz = base_xyz(xyz_f)

        self.name = name
        self.label = label

    def to_fdf(self, fout):
        fout.write("SystemName %s\n" % self.name)
        fout.write("SystemLabel %s\n" % self.label)
        fout.write("NumberOfSpecies %s\n" % self.xyz.nspecies)
        fout.write("NumberOfAtoms %s\n" % self.xyz.natom)

        cell = self.xyz.cell
        
        fout.write("%block ChemicalSpeciesLabel\n")
        for element in self.xyz.specie_labels:
            fout.write("\t%d\t%d\t%s\n" % (self.xyz.specie_labels[element], mg.Element(element).number, element))
        fout.write("%endblock ChemicalSpeciesLabel\n")
        fout.write("\n")

        fout.write("%block PAO.BasisSizes\n")
        for element in self.xyz.specie_labels:
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
        fout.write("%.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2]))
        fout.write("%.9f %.9f %.9f\n" % (cell[3], cell[4], cell[5]))
        fout.write("%.9f %.9f %.9f\n" % (cell[6], cell[7], cell[8]))
        fout.write("%endblock LatticeVectors\n")
        fout.write("\n")

        fout.write("%block AtomicCoordinatesAndAtomicSpecies\n")
        #for atom in self.xyz.atoms:
        #    fout.write("%.9f\t%.9f\t%.9f\t" % (atom.x, atom.y, atom.z))
        #    fout.write(str(self.xyz.specie_labels[atom.name]))
        #    fout.write("\n")
        for i in range(len(self.xyz.atoms)):
            fout.write("%.9f %.9f %.9f %d # %s %d\n" % (
                self.xyz.atoms[i].x, self.xyz.atoms[i].y, self.xyz.atoms[i].z,
                self.xyz.specie_labels[self.xyz.atoms[i].name],
                self.xyz.atoms[i].name,
                i+1))
            #fout.write(str(self.xyz.specie_labels[self.xyz.atoms[i].name]))
            #fout.write(" # %s %d\n" % (self.xyz.atoms[i].name, i+1))
        fout.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
        fout.write("\n")

       
