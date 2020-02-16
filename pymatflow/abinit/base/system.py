"""
in control of the system structure related parameters
"""

import pymatflow.base as base
from pymatflow.base.xyz import base_xyz


class abinit_system:
    """
    """
    def __init__(self):
        self.xyz = base_xyz()

    def to_input(self, fout):
        # fout: a file stream for writing
        cell = self.xyz.cell
        fout.write("# ==================================\n")
        fout.write("# system related setting\n")
        fout.write("# ==================================\n")
        fout.write("\n")
        fout.write("acell 1 1 1 angstrom\n\n")  # scaling with 1 means no actually scaling of rprim by acell
        fout.write("rprim\n")
        #fout.write("%.9f %.9f %.9f\n" % (cell[0], cell[1], cell[2]))
        #fout.write("%.9f %.9f %.9f\n" % (cell[3], cell[4], cell[5]))
        #fout.write("%.9f %.9f %.9f\n" % (cell[6], cell[7], cell[8]))
        for i in range(3):
            fout.write("%.9f %.9f %.9f\n" % (cell[i][0], cell[i][1], cell[i][2]))
        fout.write("\n")
        fout.write("ntypat %d\n\n" % self.xyz.nspecies)
        fout.write("natom %d\n\n" % self.xyz.natom)
        fout.write("typat\n")
        # abinit 不允许输入文件列数超过264, 因此如果原子数太多
        # 这里的typat要分多行列出
        # 利用余数设置如果一行超过30个原子就换行
        i = 0
        for atom in self.xyz.atoms:
            fout.write("%d " % self.xyz.specie_labels[atom.name])
            if i % 30 == 29:
                fout.write("\n")
            i += 1
        fout.write("\n\n")
        fout.write("znucl ")
        for element in self.xyz.specie_labels:
            fout.write(str(base.element[element].number))
            fout.write(" ")
        fout.write("\n")
        fout.write("\n")
        fout.write("xangst\n")
        for atom in self.xyz.atoms:
            fout.write("%.9f %.9f %.9f\n" % (atom.x, atom.y, atom.z))
        fout.write("\n")
        #
