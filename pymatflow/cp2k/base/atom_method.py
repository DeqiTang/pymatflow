#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

from pymatflow.cp2k.base.atom_print import cp2k_atom_print

"""
Usage:
"""



class cp2k_atom_method_external_vxc:
    def __init__(self):
        self.params = {
                }
        self.status = false


    def to_input(self, fout):
        fout.write("\t\t&EXTERNAL_VXC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END EXTERNAL_VXC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_atom_method_xc_adiabatic_rescaling:
    def __init__(self):
        self.params = {
                }
        self.status = false


    def to_input(self, fout):
        fout.write("\t\t\t\t&ADIABATIC_RESCALING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END ADIABATIC_RESCALING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_method_xc_hf_hf_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_method_xc_hf_hf_info:
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_atom_method_xc_hf_hf_info_each()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&HF_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END HF_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_atom_method_xc_hf_interaction_potential:
    def __init__(self):
        self.params = {
                }
        self.status = false

        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&INTERACTION_POTENTIAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END INTERACTION_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_atom_method_xc_hf_load_balance_print_each:
    def __init__(self):
        self.params = {
                }
        self.status = false

        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_method_xc_hf_load_balance_print:
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_atom_method_xc_hf_load_balance_print_each()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_method_xc_hf_load_balance:
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.printout = cp2k_atom_method_xc_hf_load_balance_print()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&LOAD_BALANCE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t\t&END LOAD_BALANCE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_atom_method_xc_hf_memory:
    def __init__(self):
        self.params = {
                }
        self.status = false

        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&MEMORY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END MEMORY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_atom_method_xc_hf_periodic:
    def __init__(self):
        self.params = {
                }
        self.status = false

        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&PERIODIC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END PERIODIC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_atom_method_xc_hf_screening:
    def __init__(self):
        self.params = {
                }
        self.status = false

        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&SCREENING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END SCREENING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_method_xc_hf:
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.hf_info = cp2k_atom_method_xc_hf_hf_info()
        self.interaction_potential = cp2k_atom_method_xc_hf_interaction_potential()
        self.load_balance = cp2k_atom_method_xc_hf_load_balance()
        self.memory = cp2k_atom_method_xc_hf_memory()
        self.periodic = cp2k_atom_method_xc_hf_periodic()
        self.screening = cp2k_atom_method_xc_hf_screning()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&HF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.hf_info.status == True:
            self.hf_info.to_input(fout)
        if self.interaction_potential.status == True:
            slef.interaction_potential.to_input(fout)
        if self.load_balance.status == True:
            self.load_balance.to_input(fout)
        if self.memory.status == True:
            self.memory.to_input(fotu)
        if self.periodic.status == True:
            self.periodic.to_input(fout)
        if self.screening.status == True:
            self.screening.to_input(fout)
        fout.write("\t\t\t&END HF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "HF_INFO":
                self.hf_info.set_params({item: params[item]})
            elif item.split("-")[3] == "INTERACTION_POTENTIAL":
                self.interaction_potential.set_params({item: params[item]})
            elif item.split("-")[3] == "LOAD_BALANCE":
                self.load_balance.set_params({item: params[item]})
            elif item.split("-")[3] == "MEMORY":
                self.memory.set_params({item: params[item]})
            elif item.split("-")[3] == "PERIODIC":
                self.periodic.set_params({item: params[item]})
            elif item.split("-")[3] == "SCREENING":
                self.screening.set_params({item: params[item]})
            else:
                pass



class cp2k_atom_method_xc_vdw_potential_non_local:
    def __init__(self):
        self.params = {
                }
        self.status = false


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&NON_LOCAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END NON_LOCAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_atom_method_xc_vdw_potential_pair_potential_print_dftd_each:
    def __init__(self):
        self.params = {
                }
        self.status = false


    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_method_xc_vdw_potential_pair_potential_print_dftd:
    def __init__(self):
        self.params = {
                }
        self.status = false
    
        self.each = cp2k_atom_method_xc_vdw_potential_pair_potential_print_dftd_each()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&PRINT_DFTD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END PRINT_DFTD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_atom_method_xc_vdw_potential_pair_potential:
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.print_dftd = cp2k_atom_method_xc_vdw_potential_pair_potential_print_dftd()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&PAIR_POTENTIAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.print_dftd.status == True:
            self.print_dftd.to_input(fout)
        fout.write("\t\t\t\t&END PAIR_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "PRINT_DFTD":
                self.print_dftd.set_params({item: params[item]})
            else:
                pass


class cp2k_atom_method_xc_vdw_potential:
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.non_local = cp2k_atom_method_xc_vdw_potential_non_local()
        self.pair_potential = cp2k_atom_method_xc_vdw_potential_pair_potential()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&VDW_POTENTIAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.non_local.status == True:
            self.non_local.to_input(fout)
        if self.pair_potential.status == True:
            self.pair_potential.to_input(fout)
        fout.write("\t\t\t&END VDW_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "NON_LOCAL":
                self.non_local.set_params({item: params[item]})
            elif item.split("-")[3] == "PAIR_POTENTIAL":
                self.pair_potential.set_params({item: params[item]})
            else:
                pass


class cp2k_atom_method_xc_wf_correlation:
    def __init__(self):
        self.params = {
                }
        self.status = false


    def to_input(self, fout):
        fout.write("\t\t\t\t&WF_CORRELATION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END WF_CORRELATION\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_atom_method_xc_xc_functional:
    def __init__(self):
        self.params = {
                }
        self.status = false


    def to_input(self, fout):
        fout.write("\t\t\t\t&XC_FUNCTIONAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END XC_FUNCTIONAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_atom_method_xc_xc_grid:
    def __init__(self):
        self.params = {
                }
        self.status = false


    def to_input(self, fout):
        fout.write("\t\t\t\t&XC_GRID\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END XC_GRID\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_atom_method_xc_xc_potential:
    def __init__(self):
        self.params = {
                }
        self.status = false


    def to_input(self, fout):
        fout.write("\t\t\t\t&XC_POTENTIAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END XC_POTENTIAL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_method_xc:
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.adiabatic_rescaling = cp2k_atom_method_xc_adiabatic_rescaling()
        self.hf = cp2k_atom_method_xc_hf()
        self.vdw_potential = cp2k_atom_method_xc_vdw_potential()
        self.wf_correlation = cp2k_atom_method_xc_wf_correlation()
        self.xc_functional = cp2k_atom_method_xc_xc_functional()
        self.xc_grid = cp2k_atom_method_xc_xc_grid()
        self.xc_potential = cp2k_atom_method_xc_xc_potential()
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t&XC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.adiabatic_rescaling.status == True:
            self.adiabatic_rescaling.to_input(fout)
        if self.hf.status == True:
            self.hf.to_input(fout)
        if self.vdw_potential.status == True:
            self.vde_potential.to_input(fout)
        if self.wf_correlation.status == True:
            self.wf_correlation.to_input(fout)
        if self.xc_functional.status == True:
            self.xc_functional.to_input(fout)
        if self.xc_grid.status == True:
            self.xc_grid.to_input(fout)
        if self.xc_potential.status == True:
            self.xc_potential.to_input(fout)
        fout.write("\t\t&END XC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "ADIABATIC_RESCALING":
                self.adiabatic_rescaling.set_params({item: params[item]})
            elif item.split("-")[2] == "HF":
                self.hf.set_params({item: params[item]})
            elif item.split("-")[2] == "VDW_POTENTIAL":
                self.vdw_potential.set_params({item: params[item]})
            elif item.split("-")[2] == "WF_CORRELATION":
                self.wf_correlation.set_params({item: params[item]})
            elif item.split("-")[2] == "XC_FUNCTIONAL":
                self.xc_functional.set_params({item: params[item]})
            elif item.split("-")[2] == "XC_GRID":
                self.xc_grid.set_params({item: params[item]})
            elif item.split("-")[2] == "XC_POTENTIAL":
                self.xc_potential.set_params({item: params[item]})
            else:
                pass



class cp2k_atom_method_zmp_restart:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&RESTART\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RESTART\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_atom_method_zmp:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.restart = cp2k_atom_method_zmp_restart()
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t&ZMP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.restart.status == True:
            self.restart.to_input(fout)
        fout.write("\t\t&END ZMP\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "RESTART":
                self.restart.set_params({item: params[item]})
            else:
                pass


class cp2k_atom_method:
    def __init__(self):
        self.params = {
                }
        self.status = false
        
        self.external_vxc = cp2k_atom_method_external_vxc()
        self.xc = cp2k_atom_method_xc()
        self.zmp = cp2k_atom_method_amp()
        # basic setting

    def to_input(self, fout):
        fout.write("\t&METHOD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.external_vxc.status == True:
            sefl.external_vxc.to_input(fout)
        if self.xc.status == True:
            self.xc.to_input(fout)
        if self.zmp.status == True:
            self.zmp.to_input(fout)
        fout.write("\t&END METHOD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EXTERNAL_VXC":
                self.external_vxc.set_params({item: params[item]})
            elif item.split("-")[1] == "XC":
                self.xc.set_params({item: params[item]})
            elif item.split("-")[1] == "ZMP":
                self.zmp.set_params({item: params[item]})
            else:
                pass



