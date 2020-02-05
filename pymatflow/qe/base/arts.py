#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import pymatgen as mg
import os
import sys
import re

from pymatflow.base.xyz import base_xyz


"""
usage:
"""

class qe_arts:
    """
    Control: ATOMIC_SPECIES, ATOMIC_POSITIONS, 
             K_POINTS, CELL_PARAMETERS, CONSTRAINTS, 
             OCCUPATIONS, ATOMIC_FORCES
    """
    def __init__(self):
        self.xyz = base_xyz()

        self.cell_params = {
                "cell_dynamics": None,
                "press": None,
                "wmass": None,
                "cell_factor": None,
                "press_conv_thr": None,
                "cell_dofree": None,
                }

        self.kpoints_option = "automatic"
        self.kpoints_mp = [1, 1, 1, 0, 0, 0]

        self.ifstatic = True # used to determine how to put atomic coordinates to input file
        self.atomic_forces_status = False # default is no force acting on system

    def to_in(self, fout, coordtype="angstrom"):
        # fout: a file stream for writing
        
        fout.write("ATOMIC_SPECIES\n")
        upf_all = [s for s in os.listdir("./") if s.split(".")[-1] == "UPF"]
        for element in self.xyz.specie_labels:
            for upf in upf_all:
                if upf.split(".")[0] == element:
                    pseudo_file =upf
                    break
            fout.write("%s %f %s\n" % (element, mg.Element(element).atomic_mass, pseudo_file))
            pseudo_file = None
            # after pseudo_file used, set it to None to avoid it will be used in the next element
        fout.write("\n")
        cell = self.xyz.cell
        fout.write("CELL_PARAMETERS angstrom\n")
        #fout.write("%.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2]))
        #fout.write("%.9f %.9f %.9f\n" % (cell[3], cell[4], cell[5]))
        #fout.write("%.9f %.9f %.9f\n" % (cell[6], cell[7], cell[8]))
        for i in range(3):
            fout.write("%.9f %.9f %.9f\n" % (cell[i][0], cell[i][1], cell[i][2]))
        fout.write("\n")
        if coordtype == "angstrom":
            fout.write("ATOMIC_POSITIONS angstrom\n")
            if self.ifstatic == True:
                for atom in self.xyz.atoms:
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z))
            elif self.ifstatic == False:
                for atom in self.xyz.atoms:
                    fout.write("%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z))
                    for fix in atom.fix:
                        if fix == True:
                            fout.write("\t0")
                        elif fix == False:
                            fout.write("\t1")
                    fout.write("\n")
            else:
                print("===============================================\n")
                print("warning: qe.base.arts.to_in():\n")
                print("arts.ifstatic could only be True or False\n")
                sys.exit(1)
            fout.write("\n")
        elif coordtype == "crystal":
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
            fout.write("ATOMIC_POSITIONS crystal\n")
            if self.ifstatic == True:
                for k in range(self.xyz.natom):
                    fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (self.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
            elif self.ifstatic == False:
                for k in range(self.xyz.natom):
                    fout.write("%s\t%.9f\t%.9f\t%.9f" % (self.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                    for fix in self.xyz.atoms[k].fix:
                        if fix == True:
                            fout.write("\t0")
                        elif fix == False:
                            fout.write("\t1")
                    fout.write("\n")
            else:
                print("===============================================\n")
                print("warning: qe.base.arts.to_in():\n")
                print("arts.ifstatic could only be True or False\n")
                sys.exit(1)
            fout.write("\n")
        # end crystal type ATOMIC_POSITIONS

        # writing KPOINTS to the fout
        self.write_kpoints(fout)
        # =========================
        #
        # writing forces act on atoms
        if self.atomic_forces_status == True:
            self.write_atomic_forces(fout)
        # =========================

    def write_kpoints(self, fout):
        # fout: a file stream for writing
        if self.kpoints_option == "automatic":
            fout.write("K_POINTS %s\n" % self.kpoints_option)
            fout.write("%d %d %d %d %d %d\n" % (
                self.kpoints_mp[0],
                self.kpoints_mp[1],
                self.kpoints_mp[2],
                self.kpoints_mp[3],
                self.kpoints_mp[4],
                self.kpoints_mp[5]
                ))
        elif self.kpoints_option == "gamma":
            fout.write("K_POINTS gamma\n")
        elif self.kpoints_option == "crystal_b":
            fout.write("K_POINTS %s\n" % self.kpoints_option)
            fout.write("%d\n" % len(self.crystal_b))
            for i in range(len(self.crystal_b)):
                fout.write("%f %f %f %d #%s\n" % (
                    self.crystal_b[i][0],
                    self.crystal_b[i][1],
                    self.crystal_b[i][2],
                    self.crystal_b[i][4] if self.crystal_b[i][4] != "|" else  0, # 0 is for the disconnected k point
                    self.crystal_b[i][3],
                    ))
        elif self.kpoints_option == "tpiba_b":
            pass
        

    def set_kpoints(self, kpoints_mp=[1, 1, 1, 0, 0, 0], option="automatic", crystal_b=None):
        """
        crystal_b:
            the high symmetry k point path used in bands structure calculation
            in format like this:
            
            [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]
            
            if connect_indicator in a kpoint is an integer, then it will connect to the following point
            through the number of kpoints defined by connect_indicator.
            
            if connect_indicator in a kpoint is '|', then it will not connect to the following point,
        TODO: 
        Note:
            "automatic" was controlled by kpoints_mp
            "gamma" was also handled internally
            "crystal_b" as controlled by crystal_b
        Plan:
        """
        if option == "automatic":
            self.kpoints_option = option
            self.kpoints_mp = kpoints_mp
            return
        if option == "gamma":
            self.kpoints_option = option
            return
        if option == "crystal_b":
            self.kpoints_option = option
            self.crystal_b = crystal_b
            return

    def write_atomic_forces(self, fout):
        fout.write("ATOMIC_FORCES\n")
        for i in range(len(self.xyz.atoms)):
            fout.write("%s\t%.9f\t%.9f\t%.9f\n" % (self.xyz.atoms[i].name, self.atomic_forces[i][0], self.atomic_forces[i][1] , self.atomic_forces[i][2]))

        fout.write('\n')

    def set_atomic_forces(self, pressure=None, direction=None):
        """
        set ATOMIC_FORCES
        pressure:
            in unit of Pa
        direction:
            x | y | z
        Note:
            currently only support unidirectional forces acting on all atoms of the cubic system.
            and the user provide pressure and direction of force, this function will calculate
            the corresponding force accordign to the cell.
        """

        if pressure == None or direction == None:
            self.atomic_forces_status = False
            return
        else:
            self.atomic_forces_status = True

        if direction == "x":
            area = np.sqrt(self.xyz.cell[1][0]**2 + self.xyz.cell[1][1]**2 + self.xyz.cell[1][2]**2) * np.sqrt(self.xyz.cell[2][0]**2 + self.xyz.cell[2][1]**2 + self.xyz.cell[2][2]**2) # in unit of Anstrom^2
            # 1 Hartree/Bohr = 8.238 7225(14)×10^8 N
            # 1 Ry/Bohr = 4.119358925x10^8 N
            # force is in unit of Ry/a.u.
            force = area * 1.0e-20 * pressure / (4.119358925e8)
            self.atomic_forces = np.zeros((len(self.xyz.atoms), 3))
            self.atomic_forces[:, 0] = force
    
    def basic_setting(self, ifstatic=True):
        self.ifstatic = ifstatic

    #

