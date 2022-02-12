"""
in control of the ions step related parameters
"""
import numpy as np
import sys
import os
import shutil

from pymatflow.siesta.group import SiestaVariableGroup

"""
Usage:
"""


class SiestaIons(SiestaVariableGroup):
    """
    """
    def __init__(self):
        super().__init__()
        self.incharge = [
                "WriteCoorXmol", "MD.TypeOfRun", "MD.TypeOfRun", "MD.VariableCell", "MD.ConstantVolume",
                "MD.MaxForceTol", "MD.MaxStressTol", "MD.Steps", "MD.MaxDispl", "MD.PreconditionVariableCell",
                ]

    def to_string(self):
        out = ""
        out += super().to_string()
        return out

    def basic_setting(self, option="opt"):
        """
        :param option:
            opt or md or phonon
        """
        if option == "opt":
            self.set_param("MD.TypeOfRun", "CG")   # CG, Broyden, 
            self.set_param("MD.VariableCell", "false")
            self.set_param("MD.ConstantVolume", "true")
            self.set_param("MD.MaxForceTol", 0.001, "eV/Ang")
            self.set_param("MD.MaxStressTol", 0.01, "GPa")
            self.set_param("MD.Steps", 60)
            self.set_param("MD.MaxDispl", 0.2, "Bohr")
            self.set_param("MD.PreconditionVariableCell", 5, "Ang")

            self.set_param("WriteCoorXmol", "true")
            self.set_param("WriteMDXmol", "true")

        elif option == "md":
            self.set_param("MD.TypeOfRun", "Verlet") # Verlet, Nose, ParrinelloRahman, NoseParrinelloRahman, Anneal
            self.set_param("MD.InitialTimeStep", 1)
            self.set_param("MD.FinalTimeStep", 1000) # default is MD.Steps 
            self.set_param("MD.LengthTimeStep", 1, "fs")
            self.set_param("MD.InitialTemperature", 0, "K")
            self.set_param("MD.TargetTemperature", 0, "K")

            self.params["WriteCoorXmol"] = "true"
            self.params["WriteMDXmol"] = "true"
        elif option == "phonon":
            # Phonon Calculation SETTING
            self.set_param("MD.TypeOfRun", "FC")
            self.set_param("MD.FCDispl", 0.04, "Bohr")
            self.set_param("MD.FCFirst", 1)
            self.set_param("MD.FCLast", 1)
 
        #
        #
