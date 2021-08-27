"""
in control of the system structure related parameters
"""

import numpy as np
import pymatflow.base as base
from pymatflow.base.xyz import BaseXyz


class AbinitSystem:
    """
    """
    def __init__(self):
        self.xyz = BaseXyz()
        self.status = True
        #
        self.coordtype = "cartesian" # can be cartesian or reduced

    def to_string(self, n=0):
        """
        :return input_str is the string of all the set params
        """
        cell = self.xyz.cell
        input_str = ""
        input_str += "# ==================================\n"
        input_str += "# system related setting\n"
        input_str += "# ==================================\n"
        input_str += "\n"
        input_str += "acell%s 1 1 1 angstrom\n\n" % (n if n > 0 else "")  # scaling with 1 means no actually scaling of rprim by acell
        input_str += "rprim%s\n" % (n if n > 0 else "")

        for i in range(3):
            input_str += "%.9f %.9f %.9f\n" % (cell[i][0], cell[i][1], cell[i][2])
        input_str += "\n"
        input_str += "ntypat%s %d\n\n" %(n if n > 0 else "", self.xyz.nspecies)
        input_str += "natom%s %d\n\n" % (n if n > 0 else "", self.xyz.natom)
        input_str += "typat%s\n" % (n if n > 0 else "")
        # abinit 不允许输入文件列数超过264, 因此如果原子数太多
        # 这里的typat要分多行列出
        # 利用余数设置如果一行超过30个原子就换行
        i = 0
        for atom in self.xyz.atoms:
            input_str += "%d " % self.xyz.specie_labels[atom.name]
            if i % 30 == 29:
                input_str += "\n"
            i += 1
        input_str += "\n\n"
        input_str += "znucl%s\n" % (n if n > 0 else "")
        for element in self.xyz.specie_labels:
            input_str += str(base.element[element].number)
            input_str += " "
        input_str += "\n"
        input_str += "\n"
        if self.coordtype.lower() == "cartesian":
            #input_str += "xangst%s\n" % (n if n > 0 else "")
            # xangst is disabled since v9.2 
            # now use xcart and specify the unit with angstrom
            input_str += "xcart%s\n" % (n if n > 0 else "")
            for atom in self.xyz.atoms:
                input_str += "%.9f %.9f %.9f\n" % (atom.x, atom.y, atom.z)
            input_str += "Angstrom\n"
        elif self.coordtype.lower() == "reduced":
            latcell = np.array(self.xyz.cell)
            convmat = np.linalg.inv(latcell.T)
            crystal_coord = np.zeros([self.xyz.natom, 3])
            for i in range(self.xyz.natom):
                crystal_coord[i] = convmat.dot(np.array([self.xyz.atoms[i].x, self.xyz.atoms[i].y, self.xyz.atoms[i].z]))
            #
            input_str += "xred%s\n" % (n if n > 0 else "")
            for k in range(len(crystal_coord)):
                input_str += "%.9f %.9f %.9f\n" % (crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2])
        input_str += "\n"

        return input_str
    #
