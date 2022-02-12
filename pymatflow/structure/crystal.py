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


class Crystal:
    """ an abstraction of crystal structure
    usage:
        >>> a = Crystal()
    """
    def __init__(self):
        """
        """
        self.cell = None
        self.atoms = None
        self.kpath = None

    def from_base_xyz(self, xyz):
        """
        :param basexyz: instance of pymatflow.base.xyz.BaseXyz
        """
        self.cell = xyz.cell
        self.atoms = xyz.atoms

    def from_xyz_file(self, filepath):
        """
        """
        xyz = base.BaseXyz()
        xyz.get_xyz(filepath)
        self.cell = xyz.cell
        self.atoms = xyz.atoms

    def from_cif_file(self, cif):
        """
        :param cif: filepath for cif file
        """
        import pymatflow.third.aseio
        self.cell, self.atoms = aseio.read_cif(cif)

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

    def get_fractional(self):
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
        for i in range(len(self.atoms)):
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



    def write_xyz(self, filepath):
        """
        :param filepath: output xyz file path
        """
        with open(filepath, 'w') as fout:
            fout.write("%d\n" % len(self.atoms))
            fout.write("cell: %f %f %f | %f %f %f | %f %f %f\n" % (self.cell[0][0], self.cell[0][1], self.cell[0][2], self.cell[1][0], self.cell[1][1], self.cell[1][2], self.cell[2][0], self.cell[2][1], self.cell[2][2]))
            for atom in self.atoms:
                fout.write("%s\t%f\t%f\t%f\n" % (atom.name, atom.x, atom.y, atom.z))


    def to_base_xyz(self):
        """
        :return xyz: instance of pymatflow.base.xyz.BaseXyz()
        """
        xyz = base.BaseXyz()
        xyz.file=None
        xyz.cell = self.cell
        xyz.atoms = self.atoms
        xyz.natom = len(self.atoms)
        xyz.set_species_number()
        return xyz

    def remove_atom(self, number):
        """ remove one atom from self.atoms
        :param number: an integer specifying the atom to remove
        """
        del self.atoms[number]
        self.natom = len(self.atoms)

    def remove_atoms(self, number):
        """ remove several atoms from self.atoms
        :param number: a list of integer specifying atoms to remove
            index start with 0
        """
        for i in number:
            self.atoms[i] = None
        while None in self.atoms:
            self.atoms.remove(None)
        
        self.natom = len(self.atoms)