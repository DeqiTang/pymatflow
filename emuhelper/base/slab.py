#!/usr/bin/env python
# _*_ coding: utf-8 _*_

from emuhelper.base.atom import Atom
from emuhelper.base.xyz import base_xyz

class slab_xyz(base_xyz):
    """

    """
    def __init__(self, xyz_f):
        super().__init__(xyz_f)
    
    def surf_cut(self, plane, portion):
        pass

    def put_mol(self, mol):
        """
        mol: a base_xyz object
            represeting a molecule adsorbed to the slab surface
        """
        pass
