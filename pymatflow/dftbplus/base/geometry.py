"""
    
"""
import numpy as np
import os
import sys
import re

from pymatflow.base.xyz import base_xyz

class geometry:
    def __init__(self):
        self.xyz = base_xyz()

    def to_string(self):
        strout = ""
        strout += "Geometry = {\n"
        strout += "  TypeNames = { "
        for element in self.xyz.specie_labels:
            strout += "\"%s\" "
        strout += "  }\n"
        strout += "  TypesAndCoordinates [Angstrom] = {\n"
        element_order = {}
        i = 1
        for element in self.xyz.specie_labels:
            element_order[element] = i
            i += 1
        for atom in self.xyz.atoms:
            fout.write("    %d %f %f %f\n" % (element_order[atom.name], atom.x, atom.y, atom.z))
        strout += "  }\n"
        strout += "  Periodic = Yes\n"
        strout += "  LatticeVectors [Angstrom] = {\n"
        for i in range(3):
            strout += "    %f %f %f\n" % (self.xyz.cell[i][0], self.xyz.cell[i][1], self.xyz.cell[i][2])
        strout += "  }\n"
        strout += "}\n"

