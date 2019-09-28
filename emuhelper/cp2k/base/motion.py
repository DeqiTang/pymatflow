#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

"""
Usage:
"""
# ====================
# CP2K / MOTION / BAND
# ====================
class cp2k_motion_band:
    def __init__(self):
        self.params = {
                "ALIGN_FRAMES": None,
                "BAND_TYPE": None, # CI-NEB, IT-NEB, SM
                "K_SPRING": None,
                "NPROC_REP": None,
                "NUMBER_OF_REPLICA": None,
                "POT_TYPE": None,
                "PROC_DIST_TYPE": None,
                "ROTATE_FRAMES": None,
                "USE_COLVARS": None,
                }

# ========================
# CP2K / MOTION / CELL_OPT
# ========================
class cp2k_motion_cell_opt:
    def __init__(self):
        self.params = {
                "CONSTRAINT": None,
                "EXTERNAL_PRESSURE": None,
                "KEEP_ANGLES": None,
                "KEEP_SYMMETRY": None,
                "MAX_DR": None,
                "MAX_FORCE": None,
                "MAX_ITER": None,
                "RMS_DR": None,
                "RMS_FORCE": None,
                "OPTIMIZER": None,
                "PRESSURE_TOLLERANCE": None,
                "STEP_START_VAL": None,
                "TYPE": None,
                }

class cp2k_motion_constraint:
    def __init__(self):
        self.params = {
                "CONSTRAINT_INIT": None,
                "ROLL_TOLERANCE": None,
                "SHAKE_TOLERANCE": None,
                }

# ========================
# CP2K / MOTION / DRIVER
# ========================
class cp2k_motion_driver:
    def __init__(self):
        pass

class cp2k_motion_flexible_partitioning:
    def __init__(self):
        pass

class cp2k_motion_free_energy:
    def __init__(self):
        pass

class cp2k_motion_geo_opt:
    def __init__(self):
        self.params = {
                "MAX_DR": None,
                "MAX_FORCE": None,
                "MAX_ITER": None,
                "RMS_DR": None,
                "RMS_FORCE": None,
                "OPTIMIZER": None,
                "STEP_START_VAL": None,
                "TYPE": None,
                }
 
class cp2k_motion_mc:
    def __init__(self):
        pass

class cp2k_motion_md:
    def __init__(self):
        self.params = {
                "ANGVEL_TOL": None,
                "ANGVEL_ZERO": None,
                "ANNEALING": None,
                "ANNEALING_CELL": None,
                "COMVEL_TOL": None,
                "DISPLACEMENT_TOL": None,
                "ECONS_START_VAL": None,
                "ENSEMBLE": None,
                "INITIAL_METHOD": None,
                "MAX_STEPS": None,
                "SCALE_TEMP_KIND": None,
                "STEPS": None,
                "STEP_START_VAL": None,
                "TEMPERATURE": None,
                "TEMPERATURE_ANNEALING": None,
                "TEMP_KIND": None,
                "TEMP_TOL": None,
                "TIMESTEP": None,
                "TIME_START_VAL": None,
                }

class cp2k_motion_pint:
    def __init__(self):
        self.params = {
                "DT": None,
                "FIX_CENTROID_POS": None,
                "HARM_INT": None,
                "ITERATION": None,
                "MAX_STEP": None,
                "NRESPA": None,
                "NUM_STEPS": None,
                "P": None,
                "PROC_PER_REPLICA": None,
                "PROPAGATOR": None,
                "TEMP": None,
                "TRANSFORMATION": None,
                "T_TOL": None,
                }

class cp2k_motion_print:
    def __init__(self):
        pass

class cp2k_motion_shell_opt:
    def __init__(self):
        self.params = {
                "MAX_DR": None,
                "MAX_FORCE": None,
                "MAX_ITER": None,
                "RMS_DR": None,
                "RMS_FORCE": None,
                "OPTIMIZER": None,
                "STEP_START_VAL": None,
                }

class cp2k_motion_tmc:
    def __init__(self):
        pass


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


    def to_motion(self, fname):
        with open(fname, 'a') as fout:
            fout.write("&MOTION\n")
            for item in self.params:
                fout.write("\t%s %s\n" % (item, self.params[item]))
            fout.write("&END MOTION\n")
            fout.write("\n")

