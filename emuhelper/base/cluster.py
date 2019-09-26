#!/usr/bin/env python
# _*_ coding: utf-8 _*_

from emuhelper.base.atom import Atom
from emuhelper.base.xyz import base_xyz


class cluster_xyz(base_xyz):
    """

    """
    def __init__(self, xyz_file):
        super().__init__(xyz_file)

    def build_cluster_sphere(self, radius):
        """
        radius: 球形cluster的半径, 单位是Angstrom
        """
        center = self.get_center_xyz()
        for atom in self.atoms:
            if (atom.x - center[0])**2 + (atom.y - center[1])**2 + (atom.z - center[2])**2 > radius**2:
                self.atoms.remove(atom)
        self.natom = len(self.atoms)
        
    def get_center_xyz(self):
        all_x = [atom.x for atom in self.atoms]
        all_y = [atom.y for atom in self.atoms]
        all_z = [atom.z for atom in self.atoms]
        x = (max(all_x) + min(all_x)) / 2
        y = (max(all_y) + min(all_y)) / 2
        z = (max(all_z) + min(all_z)) / 2
        #x = min(self.atoms[:].x) + (max(self.atoms[:].x) - min(self.atoms[:].x)) / 2
        #y = min(self.atoms[:].y) + (max(self.atoms[:].y) - min(self.atoms[:].y)) / 2
        #z = min(self.atoms[:].z) + (max(self.atoms[:].z) - min(self.atoms[:].z)) / 2
        print("center is: %f %f %f" % (x, y, z))
        return x, y, z
