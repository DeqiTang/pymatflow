#!/usr/bin/env pythong
# _*_ coding: utf-8 _*_

import re

class Out:
    """

    """
    def __init__(self, fname):
        self.file = fname
        self.run_type = "ENERGY_FORCE"
        self.print_level = "LOW"
        self.time_used = None
        self.eps_scf = None

    def get_input_info(self):
        with open(self.f, 'r') as fin:
            for line in fin:
                if re.match("\seps_scf", line) is not None:
                    self.eps_scf = float(line.split()[1])
                elif re.match("\sadded MOs", line) is not None:
                    self.added_mos = int(line.split()[2]) 
                elif re.match("\sElectronic temperature [K]:", line) is not None:
                    self.electronic_temperature = float(line.split()[3])
    def get_energy_info(self):
        self.get_input_info()
        if self.run_type == "ENERGY_FORCE":
            with open(self.f, 'r') as fin:
                for line in fin:
                    if re.match("\sENERGY| Total FORCE_EVAL ( QS  ) energy (a.u.):", line) is not None:
                        self.energy = float(line.split()[8])
