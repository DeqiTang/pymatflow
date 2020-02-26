
import os
import copy
import datetime
import subprocess
import matplotlib.pyplot as plt

#from pymatflow.base.atom import Atom



class dfpt_elastic_anaddb_out:
    """
    Note:
    """
    def __init__(self):
        """
        self.file:
            the output file of anaddb run for elastic constant calc.
        """
        self.file = None
        self.run_info = {}

    def get_info(self, file):
        """
        get the general information from anaddb output file
        which is now stored in self.lines
        """
        self.file = file
        with open(self.file, 'r') as fin:
            self.lines = fin.readlines()

        self.get_elastic_constant()


    def get_elastic_constant(self):
        """
        """

        self.run_info["elastic_tensor_clamped_ion"] = []
        self.run_info["elastic_tensor_relaxed_ion"] = []
        self.run_info["compliance_tensor_clamped_ion"] = []
        self.run_info["compliance_tensor_relaxed_ion"] = []

        for i in range(len(self.lines)):
            if len(self.lines[i].split()) == 0:
                continue
            if self.lines[i].split()[0] == "Elastic" and self.lines[i].split()[1] == "Tensor" and self.lines[i].split()[2] == "(clamped":
                self.run_info["elastic_tensor_clamped_ion_unit"] = self.lines[i].split("ion)")[1][:-1]
                for j in range(6):
                    self.run_info["elastic_tensor_clamped_ion"].append([
                        float(self.lines[i+j+2].split()[0]),
                        float(self.lines[i+j+2].split()[1]),
                        float(self.lines[i+j+2].split()[2]),
                        float(self.lines[i+j+2].split()[3]),
                        float(self.lines[i+j+2].split()[4]),
                        float(self.lines[i+j+2].split()[5]),
                    ])
            if self.lines[i].split()[0] == "Elastic" and self.lines[i].split()[1] == "Tensor" and self.lines[i].split()[2] == "(relaxed":
                self.run_info["elastic_tensor_relaxed_ion_unit"] = self.lines[i].split("ion)")[1][:-1]
                for j in range(6):
                    self.run_info["elastic_tensor_relaxed_ion"].append([
                        float(self.lines[i+j+3].split()[0]),
                        float(self.lines[i+j+3].split()[1]),
                        float(self.lines[i+j+3].split()[2]),
                        float(self.lines[i+j+3].split()[3]),
                        float(self.lines[i+j+3].split()[4]),
                        float(self.lines[i+j+3].split()[5]),
                    ])
            if self.lines[i].split()[0] == "Compliance" and self.lines[i].split()[1] == "Tensor" and self.lines[i].split()[2] == "(clamped":
                self.run_info["compliance_tensor_clamped_ion_unit"] = self.lines[i].split("ion)")[1][:-1]
                for j in range(6):
                    self.run_info["compliance_tensor_clamped_ion"].append([
                        float(self.lines[i+j+2].split()[0]),
                        float(self.lines[i+j+2].split()[1]),
                        float(self.lines[i+j+2].split()[2]),
                        float(self.lines[i+j+2].split()[3]),
                        float(self.lines[i+j+2].split()[4]),
                        float(self.lines[i+j+2].split()[5]),
                    ])
            if self.lines[i].split()[0] == "Compliance" and self.lines[i].split()[1] == "Tensor" and self.lines[i].split()[2] == "(relaxed":
                self.run_info["compliance_tensor_relaxed_ion_unit"] = self.lines[i].split("ion)")[1][:-1]
                for j in range(6):
                    self.run_info["compliance_tensor_relaxed_ion"].append([
                        float(self.lines[i+j+3].split()[0]),
                        float(self.lines[i+j+3].split()[1]),
                        float(self.lines[i+j+3].split()[2]),
                        float(self.lines[i+j+3].split()[3]),
                        float(self.lines[i+j+3].split()[4]),
                        float(self.lines[i+j+3].split()[5]),
                    ])
            #
