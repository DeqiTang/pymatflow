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


class octopus_system:
    """
        an abstraction of part of input block for octopus
    """
    def __init__(self):
        self.xyz = base_xyz()
        #self.pseudo = qe_pseudo()

        self.kpoints_option = "mp"
        self.kpoints_mp = [1, 1, 1, 0, 0, 0]

        self.ifstatic = True # used to determine how to put atomic coordinates to input file

    def to_inp(self, fout, coordtype="angstrom"):
        """
        :param fout: a file stream for writing
        :param coordtype: specify coordinate format, can eigher be 'angstrom' or 'reduced'
        Note:
            Note that in Octopus the origin of coordinates is in the center of the cell, 
            so the coordinates inside the cell are in the range [-0.5, 0.5).
        """
        #self.pseudo.to_in(fout=fout, xyz=self.xyz)
        fout.write("\n")
        cell = self.xyz.cell
        fout.write("\%LatticeVectors\n")
        for i in range(3):
            fout.write("%.9f | %.9f | %.9f\n" % (cell[i][0], cell[i][1], cell[i][2]))
        fout.write("\%\n")
        fout.write("\n")
        fout.write("\%LatticeParameters\n")
        fout.write("1 | 1 | 1\n")
        fout.write("\%\n")        
        fout.write("\n")
        if coordtype == "angstrom":
            fout.write("\%Coordinates\n")
            if self.ifstatic == True:
                for atom in self.xyz.atoms:
                    fout.write("\"%s\" | %.9f | %.9f | %.9f\n" % (atom.name, atom.x, atom.y, atom.z))
            elif self.ifstatic == False:
                for atom in self.xyz.atoms:
                    fout.write("\"%s\" | %.9f | %.9f | %.9f" % (atom.name, atom.x, atom.y, atom.z))
                    #for fix in atom.fix:
                    #    if fix == True:
                    #        fout.write("\t0")
                    #    elif fix == False:
                    #        fout.write("\t1")
                    if True in atom.fix:
                        fout.write(" | no")
                    else:
                        fout.write(" | yes")
                    fout.write("\n")      
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
            fout.write("\%ReducedCoordinates\n")
            if self.ifstatic == True:
                for k in range(self.xyz.natom):
                    # minus 0.5 here because  in Octopus the origin of coordinates is in the center of the cell, 
                    # so the coordinates inside the cell are in the range [-0.5, 0.5).
                    fout.write("\"%s\" | %.9f | %.9f | %.9f\n" % (self.xyz.atoms[k].name, crystal_coord[k, 0]-0.5, crystal_coord[k, 1]-0.5, crystal_coord[k, 2]-0.5))
            elif self.ifstatic == False:
                for k in range(self.xyz.natom):
                    # minus 0.5 here because  in Octopus the origin of coordinates is in the center of the cell, 
                    # so the coordinates inside the cell are in the range [-0.5, 0.5).                    
                    fout.write("\"%s\" | %.9f | %.9f | %.9f" % (self.xyz.atoms[k].name, crystal_coord[k, 0]-0.5, crystal_coord[k, 1]-0.5, crystal_coord[k, 2]-0.5))
                    #for fix in self.xyz.atoms[k].fix:
                    #    if fix == True:
                    #        fout.write("\t0")
                    #    elif fix == False:
                    #        fout.write("\t1")
                    if True in self.xyz.atoms[k].fix:
                        fout.write(" | no")
                    else:
                        fout.write(" | yes")
                    fout.write("\n")
            else:
                print("===============================================\n")
                print("warning: octopus.base.octopus_system.to_inp():\n")
                print("arts.ifstatic could only be True or False\n")
                sys.exit(1)
            fout.write("\n")
        # end reduced

        # writing KPOINTS to the fout
        self.write_kpoints(fout)
        # =========================

    def write_kpoints(self, fout):
        """
        :param fout: a file stream for writing
        """
        if self.kpoints_option == "mp":
            fout.write("\%KPointsGrid\n")
            fout.write("%d | %d | %d\n" % (
                self.kpoints_mp[0],
                self.kpoints_mp[1],
                self.kpoints_mp[2],
                ))
        fout.write("\%\n")                
        elif self.kpoints_option == "kpath":
            # there is a trick:
            # when self.kpath[i][4] == "|"
            # we set the number of k point to connect to the next high symmetry kpoint to 0
            # this is very fantastic !!!
            fout.write("\%KPointsPath\n")
            for i in range(len(self.kpath)-2):
                if self.kpath[i][4] == "|":
                    fout.write("0 | ")
                else:
                    fout.write("%d | " % self.kpath[i][4])
                fout.write("%d\n" % self.kpath[i][-2])
            #
            for i in range(len(self.kpath)):
                fout.write("%f | %f | %f #%s\n" % sefl.kpath[i][3])
            #
            #fout.write("KPointsUseSymmetries = no\n")

    def set_kpoints(self, kpoints_mp=[1, 1, 1, 0, 0, 0], option="mp", kpath=None):
        """
        :param kpath: the high symmetry k point path used in bands structure calculation
            in format like this:

            [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]

            if connect_indicator in a kpoint is an integer, then it will connect to the following point
            through the number of kpoints defined by connect_indicator.

            if connect_indicator in a kpoint is '|', then it will not connect to the following point,
        TODO:
        Note:
            "mp" means Monkhorst-Pack scheme
        """
        if option == "mp":
            self.kpoints_option = option
            self.kpoints_mp = kpoints_mp
            return
        if option == "kpath":
            self.kpoints_option = option
            self.kpath = kpath
            return

    def basic_setting(self, ifstatic=True):
        self.ifstatic = ifstatic

    #
