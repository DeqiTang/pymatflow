"""
in control of the system structure related parameters
"""
import numpy as np
import sys
import os
import shutil

import pymatflow.base as base
from pymatflow.base.atom import Atom
from pymatflow.base.xyz import BaseXyz

from pymatflow.siesta.group import SiestaVariableGroup
from . import siesta_variable_to_string

"""
Usage:
"""


class SiestaSystem(SiestaVariableGroup):
    """
    """
    def __init__(self, name="Siesta Job", label="siesta"):
        super().__init__()
        self.xyz = BaseXyz()

        self.name = name
        self.label = label

        self.incharge = [
                "PAO.FixSplitTable",
                ]        

    def to_string(self):
        out = ""
        out += "SystemName %s\n" % self.name
        out += "SystemLabel %s\n" % self.label
        out += "NumberOfSpecies %s\n" % self.xyz.nspecies
        out += "NumberOfAtoms %s\n" % self.xyz.natom

        cell = self.xyz.cell

        for item in self.params:
            if self.params[item] is not None:
                out += siesta_variable_to_string(self.params[item])
                out += "\n"

        out += "%block ChemicalSpeciesLabel\n"
        for element in self.xyz.specie_labels:
            out += "\t%d\t%d\t%s\n" % (self.xyz.specie_labels[element], base.element[element].number, element)
        out += "%endblock ChemicalSpeciesLabel\n"
        out += "\n"

        out += "%block PAO.BasisSizes\n"
        for element in self.xyz.specie_labels:
            out += "\t%s\tDZP\n" % element
        out += "%endblock PAO.BasisSizes\n"
        out += "\n"

        out += "AtomicCoordinatesFormat ScaledCartesian\n"
        # 这里可以用ScaledCartesian也可以用Ang, 因为我的LatticeConstant 设置为1Ang
        # 这样ScaledCartesian以LatticeConstant扩展后的值实际上与Ang是一样的
        out += "AtomCoorFormatOut Ang\n"
        out += "LatticeConstant 1.00000 Ang\n"
        out += "\n"

        out += "%block LatticeVectors\n"
        #out += "%.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2])
        #out += "%.9f %.9f %.9f\n" % (cell[3], cell[4], cell[5])
        #out += "%.9f %.9f %.9f\n" % (cell[6], cell[7], cell[8])
        for i in range(3):
            out += "%.9f %.9f %.9f\n" % (cell[i][0], cell[i][1], cell[i][2])
        out += "%endblock LatticeVectors\n"
        out += "\n"

        out += "%block AtomicCoordinatesAndAtomicSpecies\n"
        #for atom in self.xyz.atoms:
        #    out += "%.9f\t%.9f\t%.9f\t" % (atom.x, atom.y, atom.z)
        #    out += str(self.xyz.specie_labels[atom.name])
        #    out += "\n"
        for i in range(len(self.xyz.atoms)):
            out += "%.9f %.9f %.9f %d # %s %d\n" % (
                self.xyz.atoms[i].x, self.xyz.atoms[i].y, self.xyz.atoms[i].z,
                self.xyz.specie_labels[self.xyz.atoms[i].name],
                self.xyz.atoms[i].name,
                i+1)
            #out += str(self.xyz.specie_labels[self.xyz.atoms[i].name])
            #out += " # %s %d\n" % (self.xyz.atoms[i].name, i+1)
        out += "%endblock AtomicCoordinatesAndAtomicSpecies\n"
        out += "\n"

        return out

    def set_name_label(self, name, label):
        self.name = name
        self.label = label
    #

