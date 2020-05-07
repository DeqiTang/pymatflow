"""
Providing an abstraction of input block for octopus in control of geometric structure input
"""
import numpy as np
import os
import sys
import re

import pymatflow.base as base
from pymatflow.base.xyz import base_xyz


"""
usage:
"""

class octopus_pseudo:
    def __init__(self):
        self.dir = os.path.abspath("./") # default dir to put pseudo files

    def to_inp(self, fout, xyz):
        fout.write("\%Species\n")
        all_file = os.listdir(self.dir)
        for element in xyz.specie_labels:
            for item in all_file:
                if re.match("(%s)(.*)(psf)" % (element), item, re.IGNORECASE):
                    fout.write("\"%s\" | species_pseudo | file | \"%s\"\n" % (element, item))
                    break


class system:
    """
        an abstraction of part of input block for octopus
    """
    def __init__(self):
        self.xyz = base_xyz()
        #self.pseudo = qe_pseudo()
        
        self.ifstatic = True # used to determine how to put atomic coordinates to input file

    def to_string(self, coordtype="angstrom"):
        """
        :param fout: a file stream for writing
        :param coordtype: specify coordinate format, can eigher be 'angstrom' or 'reduced'
        Note:
            Note that in Octopus the origin of coordinates is in the center of the cell, 
            so the coordinates inside the cell are in the range [-0.5, 0.5).
        :return out -> the string
        """
        out = ""
        #self.pseudo.to_in(fout=fout, xyz=self.xyz)
        out += "\n"
        cell = self.xyz.cell
        out += "\%LatticeVectors\n"
        for i in range(3):
            out += "%.9f | %.9f | %.9f\n" % (cell[i][0], cell[i][1], cell[i][2])
        out += "\%\n"
        out += "\n"
        out += "\%LatticeParameters\n"
        out += "1 | 1 | 1\n"
        out += "\%\n"
        out += "\n"
        if coordtype == "angstrom":
            out += "\%Coordinates\n"
            if self.ifstatic == True:
                for atom in self.xyz.atoms:
                    out += "\"%s\" | %.9f | %.9f | %.9f\n" % (atom.name, atom.x, atom.y, atom.z)
            elif self.ifstatic == False:
                for atom in self.xyz.atoms:
                    out += "\"%s\" | %.9f | %.9f | %.9f" % (atom.name, atom.x, atom.y, atom.z)
                    #for fix in atom.fix:
                    #    if fix == True:
                    #        fout.write("\t0")
                    #    elif fix == False:
                    #        fout.write("\t1")
                    if True in atom.fix:
                        out += " | no"
                    else:
                        out += " | yes"
                    out += "\n"   
        elif coordtype == "reduced":
            # crystal namely fractional coordinate can be convert from cartesian coordinates
            # the conversion process is like transformation of presentation in quantum mechanics
            # the convmat is bulid to do the conversion
            #latcell = np.array(self.xyz.cell)
            #latcell = latcell.reshape(3, 3)
            latcell = np.array(self.xyz.cell)
            convmat = np.linalg.inv(latcell.T)
            crystal_coord = np.zeros([self.xyz.natom, 3])
            for i in range(self.xyz.natom):
                crystal_coord[i] = convmat.dot(np.array([self.xyz.atoms[i].x, self.xyz.atoms[i].y, self.xyz.atoms[i].z]))
            #
            out += "\%ReducedCoordinates\n"
            if self.ifstatic == True:
                for k in range(self.xyz.natom):
                    # minus 0.5 here because  in Octopus the origin of coordinates is in the center of the cell, 
                    # so the coordinates inside the cell are in the range [-0.5, 0.5).
                    out += "\"%s\" | %.9f | %.9f | %.9f\n" % (self.xyz.atoms[k].name, crystal_coord[k, 0]-0.5, crystal_coord[k, 1]-0.5, crystal_coord[k, 2]-0.5)
            elif self.ifstatic == False:
                for k in range(self.xyz.natom):
                    # minus 0.5 here because  in Octopus the origin of coordinates is in the center of the cell, 
                    # so the coordinates inside the cell are in the range [-0.5, 0.5).                    
                    out += "\"%s\" | %.9f | %.9f | %.9f" % (self.xyz.atoms[k].name, crystal_coord[k, 0]-0.5, crystal_coord[k, 1]-0.5, crystal_coord[k, 2]-0.5)
                    #for fix in self.xyz.atoms[k].fix:
                    #    if fix == True:
                    #        fout.write("\t0")
                    #    elif fix == False:
                    #        fout.write("\t1")
                    if True in self.xyz.atoms[k].fix:
                        out += " | no"
                    else:
                        out += " | yes"
                    out += "\n"
            else:
                print("===============================================\n")
                print("warning: octopus.base.octopus_system.to_inp():\n")
                print("arts.ifstatic could only be True or False\n")
                sys.exit(1)
            out += "\n"
        # end reduced
        return out

    def basic_setting(self, ifstatic=True):
        self.ifstatic = ifstatic

    #
