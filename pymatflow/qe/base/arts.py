"""
Providing an abstraction of input block for pwscf in control of ATOMIC_SPECIES, ATOMIC_POSITIONS, KPOINTS,
CELL_PARAMETERS, CONSTRAINTS, OCCUPATIONS, ATOMIC_FORCES
"""
import numpy as np
import os
import sys
import re

import pymatflow.base as base
from pymatflow.base.xyz import BaseXyz

from pymatflow.qe.group import QeVariableGroup
from . import qe_variable_to_string

"""
usage:
"""

class QePseudo:
    def __init__(self):
        self.dir = os.path.abspath("./") # default dir to put pseudo files

    def to_in(self, fout, xyz):
        fout.write(self.to_string(xyz))

    def to_string(self, xyz):
        out = ""
        out += "ATOMIC_SPECIES\n"
        all_file = os.listdir(self.dir)
        for element in xyz.specie_labels:
            for item in all_file:
                if re.match("(%s)(.*)(upf)" % (element), item, re.IGNORECASE):
                    out += "%s %f %s\n" % (element, base.element[element].mass, item)
                    break
        return out

class QeArts(QeVariableGroup):
    """
        an abstraction of part of input block for pwscf
    """
    def __init__(self):
        super().__init__()
        self.xyz = BaseXyz()
        self.pseudo = QePseudo()

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
        """
        :param fout: a file stream for writing
        :param coordtype: specify coordinate format, can eigher be 'angstrom' or 'crystal'
        """
        fout.write(self.to_string(coordtype=coordtype))

    def to_string(self, coordtype="angstrom"):
        out = ""
        out += self.pseudo.to_string(xyz=self.xyz)
        out += "\n"
        cell = self.xyz.cell
        out += "CELL_PARAMETERS angstrom\n"
        #fout.write("%.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2]))
        #fout.write("%.9f %.9f %.9f\n" % (cell[3], cell[4], cell[5]))
        #fout.write("%.9f %.9f %.9f\n" % (cell[6], cell[7], cell[8]))
        for i in range(3):
            out += "%.9f %.9f %.9f\n" % (cell[i][0], cell[i][1], cell[i][2])
        out += "\n"
        if coordtype == "angstrom":
            out += "ATOMIC_POSITIONS angstrom\n"
            if self.ifstatic == True:
                for atom in self.xyz.atoms:
                    out += "%s\t%.9f\t%.9f\t%.9f\n" % (atom.name, atom.x, atom.y, atom.z)
            elif self.ifstatic == False:
                for atom in self.xyz.atoms:
                    out += "%s\t%.9f\t%.9f\t%.9f" % (atom.name, atom.x, atom.y, atom.z)
                    for fix in atom.fix:
                        if fix == True:
                            out += "\t0"
                        elif fix == False:
                            out += "\t1"
                    out += "\n"
            else:
                print("===============================================\n")
                print("warning: qe.base.arts.to_in():\n")
                print("arts.ifstatic could only be True or False\n")
                sys.exit(1)
            out += "\n"
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
            out += "ATOMIC_POSITIONS crystal\n"
            if self.ifstatic == True:
                for k in range(self.xyz.natom):
                    out += "%s\t%.9f\t%.9f\t%.9f\n" % (self.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2])
            elif self.ifstatic == False:
                for k in range(self.xyz.natom):
                    out += "%s\t%.9f\t%.9f\t%.9f" % (self.xyz.atoms[k].name, crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2])
                    for fix in self.xyz.atoms[k].fix:
                        if fix == True:
                            out += "\t0"
                        elif fix == False:
                            out += "\t1"
                    out += "\n"
            else:
                print("===============================================\n")
                print("warning: qe.base.arts.to_in():\n")
                print("arts.ifstatic could only be True or False\n")
                sys.exit(1)
            out += "\n"
        # end crystal type ATOMIC_POSITIONS

        # writing KPOINTS to the fout
        out += self.to_string_kpoints()
        # =========================
        #
        # writing forces act on atoms
        if self.atomic_forces_status == True:
            out += self.to_string_atomic_forces()
        # =========================
        return out

    def write_kpoints(self, fout):
        """
        :param fout: a file stream for writing
        """
        fout.write(self.to_string_kpionts())

    def to_string_kpoints(self):
        out = ""
        if self.kpoints_option == "automatic":
            out += "K_POINTS %s\n" % self.kpoints_option
            out += "%d %d %d %d %d %d\n" % (
                self.kpoints_mp[0],
                self.kpoints_mp[1],
                self.kpoints_mp[2],
                self.kpoints_mp[3],
                self.kpoints_mp[4],
                self.kpoints_mp[5]
                )
        elif self.kpoints_option == "gamma":
            out += "K_POINTS gamma\n"
        elif self.kpoints_option == "crystal_b":
            # there is a trick:
            # when self.crystal_b[i][4] == "|"
            # we set the number of k point to connect to the next high symmetry kpoint to 0
            # then after the pw.x calculation and bands.x calculation, we can see from the
            # output of bands.x, the two nieghbor high symmetry kpoints
            # have the same x coordinates, and that can be processed by qe.post.bands and
            # the corresponding post-qe-bands.py, and the label in the processed band strucutre
            # image would be in format of 'label|label', for example 'K|U'
            # this is very fantastic !!!
            out += "K_POINTS %s\n" % self.kpoints_option
            out += "%d\n" % len(self.crystal_b)
            for i in range(len(self.crystal_b)):
                out += "%f %f %f %d #%s\n" % (
                    self.crystal_b[i][0],
                    self.crystal_b[i][1],
                    self.crystal_b[i][2],
                    self.crystal_b[i][4] if self.crystal_b[i][4] != "|" else  0, # 0 is for the disconnected k point
                    self.crystal_b[i][3],
                    )
        elif self.kpoints_option == "tpiba_b":
            pass
        
        return out

    def set_kpoints(self, kpoints_mp=[1, 1, 1, 0, 0, 0], option="automatic", crystal_b=None):
        """
        :param crystal_b: the high symmetry k point path used in bands structure calculation
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
        fout.write(self.to_string_atomic_forces)

    def to_string_atomic_forces(self):
        out = ""
        out += "ATOMIC_FORCES\n"
        for i in range(len(self.xyz.atoms)):
            out += "%s\t%.9f\t%.9f\t%.9f\n" % (self.xyz.atoms[i].name, self.atomic_forces[i][0], self.atomic_forces[i][1] , self.atomic_forces[i][2])

        out += '\n'
        return out

    def set_atomic_forces(self, pressure=None, direction=None):
        """set ATOMIC_FORCES
        :param pressure:
            in unit of Pa
        :param direction:
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
