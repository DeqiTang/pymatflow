#!/usr/bin/env python
# _*_ coding: utf-8


class VaspPotcar:
    """
    """
    def __init__(self, poscar):
        self.elements = [i for i in poscar.xyz.specie_labels]
        self.potdir = ""

    def to_potcar(self, fname="POTCAR"):
        pass
