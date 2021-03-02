

import numpy as np
import sys
import os
import shutil
import pymatflow.base as base

from pymatflow.base.xyz import BaseXyz

"""
Usage:
"""

class cp2k_subsys_cell:
    def __init__(self):
        self.params = {
                "ALPHA_BETA_GAMMA": None,
                "CELL_FILE_FORMAT": None,
                "CELL_FILE_NAME": None,
                "MULTIPLE_UNIT_CELL": None,
                "PERIODIC": None,
                "SYMMETRY": None,
                }
        self.status = False

    def to_input(self, fout, xyz):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&CELL\n")
        fout.write("\t\t\tA %.9f %.9f %.9f\n" % (xyz.cell[0][0], xyz.cell[0][1], xyz.cell[0][2]))
        fout.write("\t\t\tB %.9f %.9f %.9f\n" % (xyz.cell[1][0], xyz.cell[1][1], xyz.cell[1][2]))
        fout.write("\t\t\tC %.9f %.9f %.9f\n" % (xyz.cell[2][0], xyz.cell[2][1], xyz.cell[2][2]))
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END CELL\n")
    
    def set_params(self, params):
        """
        """
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_subsys_colvar:
    def __init__(self):
        self.status = False

class cp2k_subsys_coord:
    def __init__(self):
        self.params = {
            "DEFAULT_KEYWORD": None,
            "SCALED": None,
            "UNIT": None,
        }
        self.status = False

    def to_input(self, fout, xyz):
        fout.write("\t\t&COORD\n")
        for atom in xyz.atoms:
            fout.write("\t\t\t%s %.9f %.9f %.9f\n" % (atom.name, atom.x, atom.y, atom.z))
        fout.write("\t\t&END COORD\n")

class cp2k_subsys_core_coord:
    def __init__(self):
        self.status = False

class cp2k_subsys_core_velocity:
    def __init__(self):
        self.status = False

class cp2k_subsys_kind:
    def __init__(self):
        self.status = False

        self.basis_set = dict()
        self.potential = dict()
        for i in base.element:
            self.basis_set[str(i)] = 'DZVP-MOLOPT-SR-GTH'
            self.potential[str(i)] = 'GTH-PBE'

    def to_input(self, fout, xyz):
        """
        fout: a file stream for writing
        """
        for element in xyz.specie_labels:
                fout.write("\t\t&KIND %s\n" % element)
                fout.write("\t\t\tBASIS_SET %s\n" % self.basis_set[element])
                fout.write("\t\t\tPOTENTIAL %s\n" % self.potential[element])
                fout.write("\t\t&END KIND\n")

class cp2k_subsys_multipoles:
    def __init__(self):
        self.status = False

class cp2k_subsys_print:
    def __init__(self):
        self.status = False

class cp2k_subsys_rng_init:
    def __init__(self):
        self.status = False

class cp2k_subsys_shell_coord:
    def __init__(self):
        self.status = False

class cp2k_subsys_shell_velocity:
    def __init__(self):
        self.status = False

class cp2k_subsys_topology:
    def __init__(self):
        self.status = False

    def to_input(self, fout, xyz):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&TOPOLOGY\n")
        fout.write("\t\t\tCOORD_FILE_FORMAT xyz\n")
        #fout.write("\t\t\tCOORD_FILE_NAME %s\n" % xyz.file)
        fout.write("\t\t\tCOORD_FILE_NAME %s\n" % os.path.basename(xyz.file))
        fout.write("\t\t&END TOPOLOGY\n")

class cp2k_subsys_velocity:
    def __init__(self):
        self.status = False




# ===========================
# CP2K / FORCE_EVAL / SUBSYS
# ===========================
class cp2k_subsys:
    """

    """
    def __init__(self):
        #super().__init__()
        self.xyz = BaseXyz()
        self.params = {
                "SEED": None,
                }
        self.status = False

        self.cell = cp2k_subsys_cell()

        self.colvar = cp2k_subsys_colvar()

        self.coord = cp2k_subsys_coord()

        self.core_coord = cp2k_subsys_core_coord()

        self.core_velocity = cp2k_subsys_core_velocity()

        self.kind = cp2k_subsys_kind()

        self.multipoles = cp2k_subsys_multipoles()

        self.printout = cp2k_subsys_print()

        self.rng_init = cp2k_subsys_rng_init()

        self.shell_coord = cp2k_subsys_shell_coord()

        self.shell_velocity = cp2k_subsys_shell_velocity()

        self.topology = cp2k_subsys_topology()

        self.velocity = cp2k_subsys_velocity()

        # basic_setting
        self.kind.status = True
        self.cell.status = True
        self.topology.status = True

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&SUBSYS\n")
        if self.kind.status == True:
            self.kind.to_input(fout, self.xyz)
        if self.cell.status == True:
            self.cell.to_input(fout, self.xyz)
        if self.coord.status == True:
            self.coord.to_input(fout, self.xyz)
        if self.topology.status == True:
            self.topology.to_input(fout, self.xyz)
        fout.write("\t&END SUBSYS\n")

    def set_params(self, params):
        """
            Note: parameters for sub section(like md), are handled over to sub section controllers.
        """
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "CELL":
                self.cell.set_params({item: params[item]})
            else:
                pass