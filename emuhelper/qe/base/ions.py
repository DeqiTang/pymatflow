#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import pymatgen as mg

"""
usage:
"""

class qe_ions:
    """

    """
    def __init__(self):
        self.params = {
                "ion_dynamics": None,
                "ion_positions": None,
                "pot_extrapolation": None,
                "wfc_extrapolation": None,
                "remove_rigid_rot": None,
                "ion_temperature": None,
                "tempw": None,
                "tolp": None,
                "delta_t": None,
                "nraise": None,
                "refold_pos": None,
                "upscale": None,
                "bfgs_ndim": None,
                "trust_radius_max": None,
                "trust_radius_min": None,
                "trust_radius_ini": None,
                "w_1": None,
                "w_2": None,
                }
    def to_in(self, fout):
        # fout: a file stream for writing
        fout.write("&ions\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item, str(self.params[item])))
        fout.write("/\n")
        fout.write("\n")
