"""
module for crystal structure
"""

import numpy as np
import sys
import os
import shutil
import copy
import pymatflow.base as base
from pymatflow.base.atom import Atom

"""
Usage:
"""


class crystal():
    """ an abstraction of crystal structure
    usage:
        >>> a = crystal()
    """
    def __init__(self):
        """
        Note: use 'have a' rather than 'is kind of' to make crystal() more independent
            of base_xyz().
        """
        self.cell = None
        self.atoms = None

    def from_base_xyz(self, basexyz):
        """
        :param basexyz: instance of pymatflow.base.xyz.base_xyz
        """
        self.cell = basexyz.cell
        self.atoms = basexyz.atoms

    def from_xyz_file(self, xyz):
        """
        """
        pass

    def from_cif_file(self, cif):
        """
        """
        pass

    def get_cell(self, cell):
        """
        :params cell: [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]] in unit of Anstrom
        """
        self.cell = cell

    def get_atoms(self, atoms):
        """
        :params cell: [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]] in unit of Anstrom
        :params atoms (in cartesian coordinates and unit of Anstrom)
                [
                    ["C", 0.00000, 0.0000000, 0.0000],
                    ["O", 1.12300, 3.3250000, 2.4893],
                    ....
                ]
        """
        self.atoms = [Atom(atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3]) for i in range(len(atoms))]

    def get_cell_atoms(self, cell, atoms):
        """
        :params cell: [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]] in unit of Anstrom
        :params atoms (in cartesian coordinates and unit of Anstrom)
                [
                    ["C", 0.00000, 0.0000000, 0.0000],
                    ["O", 1.12300, 3.3250000, 2.4893],
                    ....
                ]
        """
        self.cell = cell
        self.atoms = [Atom(atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3]) for i in range(len(atoms))]

    def cell(self):
        """
        :return cell: cell parameters of the structure
            [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]]
        """
        return self.cell

    def cartesian(self):
        """
        :return cartesian coordinates
            [
                ["C", 0.00000, 0.0000000, 0.0000],
                ["O", 1.12300, 3.3250000, 2.4893],
                ....
            ]

        """
        return [[self.atoms[i].name, self.atoms[i].x, self.atoms[i].y, self.atoms[i].z] for i in range(self.natom)]

    def fractional(self):
        """
        :return out: fractional coordinates
            [
                ["C", 0.00000, 0.000000, 0.0000],
                ["O", 0.00000, 0.500000, 0.0000],
                ....
            ]

        """
        out = []
        latcell = np.array(self.cell)
        convmat = np.linalg.inv(latcell.T)
        for i in range(self.natom):
            atom = []
            atom.append(self.atoms[i].name)
            atom = atom + list(convmat.dot(np.array([self.atoms[i].x, self.atoms[i].y, self.atoms[i].z])))
            out.append(atom)
        #
        return out

    def volume(self):
        """
        :return volume in unit of Angstrom^3
        """
        return np.linalg.det(self.cell)

    def build_supercell(self, n):
        """
        :param n: [n1, n2, n3]
        :return out:
            {
                "cell": [[], [], []],
                "atoms": [
                        ["C", 0.00000, 0.000000, 0.0000],
                        ["O", 0.00000, 0.500000, 0.0000],
                        ...
                    ]
            }
        Note: will not affect status of self
        """
        #
        cell = copy.deepcopy(self.cell)
        for i in range(3):
            for j in range(3):
                cell[i][j] = n[i] * self.cell[i][j]
        atoms = copy.deepcopy(self.atoms)
        # build supercell: replica in three vector one by one
        for i in range(3):
            natom_now = len(atoms)
            for j in range(n[i] - 1):
                for atom in atoms[:natom_now]:
                    x = atom.x + float(j + 1) * self.cell[i][0]
                    y = atom.y + float(j + 1) * self.cell[i][1]
                    z = atom.z + float(j + 1) * self.cell[i][2]
                    atoms.append(Atom(atom.name, x, y, z))
        return {"cell": cell, "atoms": [[atom.name, atom.x, atom.y, atom.z] for atom in atoms]}

    def to_base_xyz(self):
        """
        :return xyz: instance of pymatflow.base.xyz.base_xyz()
        """
        xyz = base.base_xyz()
        xyz.file=None
        xyz.cell = self.cell
        xyz.atoms = self.atoms
        xyz.natom = len(self.atoms)
        xyz.set_species_number()
        return xyz
