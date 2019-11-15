#!/usr/bin/env python
# _*_ coding: utf-8 _*_

from emuhelper.base.atom import Atom
from emuhelper.base.xyz import base_xyz

class slab:
    """
    
    """
    def __init__(self, xyz=None):
        """
        xyz:  an instance of base_xyz
        """
        if xyz is not None:
            self.xyz = xyz
    
    def return_xyz(self):
        """
        return the final constructed structure(an instance of base_xyz)
        """
        return self.xyzdone

    def get_xyz(self, xyz):
        """
        xyz: an instance of base_xyz
        """
        self.xyz = xyz

    def surf_cut(self, plane, portion):
        pass

    def put_mol(self, mol):
        """
        mol: a base_xyz object
            represeting a molecule adsorbed to the slab surface
        """
        pass
