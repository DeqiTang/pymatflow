#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil



from pymatflow.vasp.base.xyz import base_xyz

"""
Usage:

    make sure the xyz structure file and the POTCAR is in the directory.
    make sure the element in the POTCAR is in order of increasing the atom number of
    the element: 确保POTCAR中元素排列的顺序是原子序数从小到大
Note:

"""

class vasp_poscar:
    """
    """
    def __init__(self):
        self.xyz = base_xyz()
        self.selective_dynamics = False

    def to_poscar(self, fout, coordtype="Cartesian"):
        """
        :param fout: a file stream for write
        :param coordtype:
            Direct, DIRECT, direct, or
            Cartesian, CARTESIAN, cartesian
        """
        cell = self.xyz.cell
        fout.write("general comment\n")
        fout.write("1.0\n") # universal scaling parameters
        for vec in cell:
            fout.write("%f %f %f\n" % (vec[0], vec[1], vec[2])) # lattice vector a(3)
        for element in self.xyz.specie_labels:
            fout.write("%s " % element)
        fout.write("\n")
        for element in self.xyz.specie_labels:
            n_atom_tmp = 0
            for atom in self.xyz.atoms:
                if atom.name == element:
                    n_atom_tmp += 1
            fout.write("%d " % n_atom_tmp)
        fout.write("\n")
        # note: the order of atoms in the POSCAR might be different from
        # the input structure file xxx.xyz
        # as in vasp, we have to input it according to the element type.
        if self.selective_dynamics == True:
            fout.write("Selective dynamics\n")
            if coordtype.lower() == "cartesian":
                fout.write("Cartesian\n")
                for element in self.xyz.specie_labels:
                    for atom in self.xyz.atoms:
                        if atom.name == element:
                            fout.write("%f %f %f" % (atom.x, atom.y, atom.z))
                            for fix in atom.fix:
                                if fix == True:
                                    fout.write("\tF")
                                elif fix == False:
                                    fout.write("\tT")
                            fout.write("\n")
            elif coordtype.lower() == "direct":
                # crystal namely fractional coordinate can be convert from cartesian coordinates
                # the conversion process is like transformation of presentation in quantum mechanics
                # the convmat is bulid to do the conversion
                latcell = np.array(self.xyz.cell)
                latcell = latcell.reshape(3, 3)
                convmat = np.linalg.inv(latcell.T)
                crystal_coord = np.zeros([self.xyz.natom, 3])
                for i in range(self.xyz.natom):
                    crystal_coord[i] = convmat.dot(np.array([self.xyz.atoms[i].x, self.xyz.atoms[i].y, self.xyz.atoms[i].z]))
                #
                fout.write("Direct\n")
                for element in self.xyz.specie_labels:
                    for k in range(self.xyz.natom):
                        if self.xyz.atoms[k].name == element:
                            fout.write("%.9f\t%.9f\t%.9f" % (crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))
                            for fix in self.xyz.atoms[k].fix:
                                if fix == True:
                                    fout.write("\tF")
                                elif fix == False:
                                    fout.write("\tT")
                            fout.write("\n")
        elif self.selective_dynamics == False:
            if coordtype.lower() == "cartesian":
                fout.write("Cartesian\n")
                for element in self.xyz.specie_labels:
                    for atom in self.xyz.atoms:
                        if atom.name == element:
                            fout.write("%f %f %f\n" % (atom.x, atom.y, atom.z))
            elif coordtype.lower() == "direct":
                fout.write("Direct\n")
                # crystal namely fractional coordinate can be convert from cartesian coordinates
                # the conversion process is like transformation of presentation in quantum mechanics
                # the convmat is bulid to do the conversion
                latcell = np.array(self.xyz.cell)
                latcell = latcell.reshape(3, 3)
                convmat = np.linalg.inv(latcell.T)
                crystal_coord = np.zeros([self.xyz.natom, 3])
                for i in range(self.xyz.natom):
                    crystal_coord[i] = convmat.dot(np.array([self.xyz.atoms[i].x, self.xyz.atoms[i].y, self.xyz.atoms[i].z]))
                #
                for element in self.xyz.specie_labels:
                    for k in range(self.xyz.natom):
                        if self.xyz.atoms[k].name == element:
                            fout.write("%.9f\t%.9f\t%.9f\n" % (crystal_coord[k, 0], crystal_coord[k, 1], crystal_coord[k, 2]))

        else:
            print("===============================================\n")
            print("warning: vasp.base.xyz.vasp_xyz.to_poscar():\n")
            print("vasp.base.poscar.ifstatic could only be True or False\n")
            sys.exit(1)
        fout.write("\n")
        # end for the coordinate of the atoms

    def get_cell(self):
        """
        cell defined in xxx.xyz must be in format like this:
        cell: 4.08376 0.00000 0.00000 | 0.00000 4.00251 0.00000 | -0.05485 0.00000 8.16247
        """
        with open(self.file, 'r') as fin:
            fin.readline()
            line = fin.readline()
        return [float(line.split()[i]) for i in [1, 2, 3, 5, 6, 7, 9, 10, 11]]

    def update(self, newxyzfile):
        self.file = newxyzfile
        self.natom = 0
        self.nspecies = 0
        self.atoms = []
        self.specie_labels = dict()
        self.get_info()
        self.cell = self.get_cell()
