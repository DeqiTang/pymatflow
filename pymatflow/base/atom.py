"""
    a representation of atom with xyz coordinates
"""

import numpy as np
import sys
import os
import shutil


"""
Usage:
"""


class Atom:
    """
        a representation of atom with xyz coordinates
    """
    def __init__(self, name=None, x=0, y=0, z=0):
        self.name = name
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

        self.fix = [False, False, False] # whether fix the atom when doing geo opt or md

    def set_name(self, name):
        self.name = name
    def set_x(self, x):
        self.x = float(x)
    def set_y(self, y):
        self.y = float(y)
    def set_z(self, z):
        self.z = float(z)

