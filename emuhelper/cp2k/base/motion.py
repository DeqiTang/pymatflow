#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from emuhelper.cp2k.base.motion_print import cp2k_motion_print
from emuhelper.cp2k.base.motion_band import cp2k_motion_band
from emuhelper.cp2k.base.motion_cell_opt import cp2k_motion_cell_opt
from emuhelper.cp2k.base.motion_geo_opt import cp2k_motion_geo_opt
from emuhelper.cp2k.base.motion_constraint import cp2k_motion_constraint
from emuhelper.cp2k.base.motion_driver import cp2k_motion_driver
from emuhelper.cp2k.base.motion_free_energy import cp2k_motion_free_energy
from emuhelper.cp2k.base.motion_mc import cp2k_motion_mc
from emuhelper.cp2k.base.motion_md import cp2k_motion_md
from emuhelper.cp2k.base.motion_pint import cp2k_motion_pint
from emuhelper.cp2k.base.motion_shell_opt import cp2k_motion_shell_opt
from emuhelper.cp2k.base.motion_tmc import cp2k_motion_tmc
from emuhelper.cp2k.base.motion_flexible_partitioning import cp2k_motion_flexible_partitioning


"""
Usage:
"""

class cp2k_motion:
    """
    """
    def __init__(self, project_name="Ab-initio", run_type="NONE"):
        self.params = {
                }

        self.band = cp2k_motion_band()

        self.cell_opt = cp2k_motion_cell_opt()
        
        self.constraint = cp2k_motion_constraint() 
        
        self.driver = cp2k_motion_driver()

        self.flexible_partitioning = cp2k_motion_flexible_partitioning()

        self.free_energy = cp2k_motion_free_energy()

        self.geo_opt = cp2k_motion_geo_opt()
       
        self.mc = cp2k_motion_mc()

        self.md = cp2k_motion_md()
        
        self.pint = cp2k_motion_pint()
        
        self.printout = cp2k_motion_print() 
                
        self.shell_opt = cp2k_motion_shell_opt()
        
        self.tmc = cp2k_motion_tmc() 

        self.run_type = "GEO_OPT"

    def to_input(self, fout):
        fout.write("&MOTION\n")
        for item in self.params:
            fout.write("\t%s %s\n" % (item, self.params[item]))
        if self.run_type == "GEO_OPT":
            self.geo_opt.to_motion(fout)
        elif self.run_type == "CELL_OPT":
            self.cell_opt.to_motion(fout)
        elif self.run_type == "MC":
            self.mc.to_motion(fout)
        elif self.run_type == "MD":
            self.md.to_motion(fout)
        elif self.run_type == "PINT":
            self.pint.to_motion(fout)
        elif self.run_type == "SHELL_OPT":
            self.shell_opt.to_motion(fout)
        elif self.run_type == "BAND":
            self.band.to_motion(fout)

        self.printout.to_motion(fout)

        fout.write("&END MOTION\n")
        fout.write("\n")
    
    def set_type(self, run_type):
        """
        runtype: CELL_OPT, GEO_OPT, MC, MD, SHELL_OPT
        """
        self.run_type = run_type

    def set_params(self, params):
        """
        parameters for sub section(like md), are handled over 
        to sub section controllers.
        """
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "MD":
                self.md.set_params({item: params[item]})
            elif item.split("-")[0] == "BAND":
                self.band.set_params({item: params[item]})
            elif item.split("-")[0] == "CELL_OPT":
                self.cell_opt.set_params({item: params[item]})
            elif item.split("-")[0] == "GEO_OPT":
                self.geo_opt.set_params({item: params[item]})
            elif item.split("-")[0] == "MC":
                self.mc.set_params({item: params[item]})
            elif item.split("-")[0] == "TMC":
                self.tmc.set_params({item: params[item]})
            elif item.split("-")[0] == "SHELL_OPT":
                self.shell_opt.set_params({item: params[item]})
            elif item.split("-")[0] == "CONSTRAINT":
                self.constraint.set_params({item: params[item]})
            elif item.split("-")[0] == "DRIVER":
                self.driver.set_params({item: params[item]})
            elif item.split("-")[0] == "PINT":
                self.pint.set_params({item: params[item]})
            elif item.split("-")[0] == "FLEXIBLE_PARTITIONING":
                self.flexible_partitioning.set_params({item: params[item]})
            elif item.split("-")[0] == "FREE_ENERGY":
                self.free_energy.set_params({item: params[item]})
            elif item.split("-")[0] == "PRINT":
                self.printout.set_params({item: params[item]})
