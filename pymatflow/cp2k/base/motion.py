"""
a representation for MOTION
"""
import numpy as np
import sys
import os
import shutil

from pymatflow.cp2k.base.motion_print import cp2k_motion_print
from pymatflow.cp2k.base.motion_band import cp2k_motion_band
from pymatflow.cp2k.base.motion_cell_opt import cp2k_motion_cell_opt
from pymatflow.cp2k.base.motion_geo_opt import cp2k_motion_geo_opt
from pymatflow.cp2k.base.motion_constraint import cp2k_motion_constraint
from pymatflow.cp2k.base.motion_driver import cp2k_motion_driver
from pymatflow.cp2k.base.motion_free_energy import cp2k_motion_free_energy
from pymatflow.cp2k.base.motion_mc import cp2k_motion_mc
from pymatflow.cp2k.base.motion_md import cp2k_motion_md
from pymatflow.cp2k.base.motion_pint import cp2k_motion_pint
from pymatflow.cp2k.base.motion_shell_opt import cp2k_motion_shell_opt
from pymatflow.cp2k.base.motion_tmc import cp2k_motion_tmc
from pymatflow.cp2k.base.motion_flexible_partitioning import cp2k_motion_flexible_partitioning


"""
Usage:
"""

class cp2k_motion:
    """
    """
    def __init__(self, project_name="ab-initio", run_type="NONE"):
        self.params = {
                }

        self.status = False

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

        # basic setting
        self.printout.status = True

    def to_input(self, fout):
        fout.write("&MOTION\n")
        for item in self.params:
            fout.write("\t%s %s\n" % (item, self.params[item]))
        if  self.geo_opt.status == True:
            self.geo_opt.to_input(fout)
        if  self.cell_opt.status == True:
            self.cell_opt.to_input(fout)
        if  self.mc.status == True:
            self.mc.to_input(fout)
        if self.md.status == True:
            self.md.to_input(fout)
        if self.pint.status == True:
            self.pint.to_input(fout)
        if self.shell_opt.status == True:
            self.shell_opt.to_input(fout)
        if self.band.status == True:
            self.band.to_input(fout)
        if self.free_energy.status == True:
            self.free_energy.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("&END MOTION\n")
        fout.write("\n")
    
    def set_type(self, run_type):
        """
        :param runtype: CELL_OPT, GEO_OPT, MC, MD, SHELL_OPT
        """
        if run_type == "GEO_OPT":
            self.geo_opt.status = True
        if run_type == "CELL_OPT":
            self.cell_opt.status = True
            self.check_cell_opt()
        if run_type == "MC":
            self.mc.status = True
        if run_type == "MD":
            self.md.status = True
        if run_type == "PINT":
            self.pint.status = True
        if run_type == "SHELL_OPT":
            self.shell_opt.status = True
        if run_type == "BAND":
            self.band.status = True

    def set_params(self, params):
        """
            Note: parameters for sub section(like md), are handled over to sub section controllers.
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
            else:
                pass
        # need to do some checking
        self.check_cell_opt()

    def check_cell_opt(self):
        """
            Note: whenever we set cell opt params whether through self.set_type or self.set_params
                we should check_cell_opt() to make sure some rules is followed. see below
        """
        # if the MOTION-CELL_OPT-TYPE == "GEO_OPT"
        # we also need to provide GEO_OPT section
        if self.cell_opt.params["TYPE"] is not None and self.cell_opt.params["TYPE"].upper() == "GEO_OPT":
            self.geo_opt.status ==  True

