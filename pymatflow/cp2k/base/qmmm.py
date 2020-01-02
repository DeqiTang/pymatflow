#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg


"""
usage:
"""

# ============================================
# CP2K / QMMM
#=============================================

class cp2k_qmmm_cell_cell_ref:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CELL_REF\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END CELL_REF\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_cell:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cell_ref = cp2k_qmmm_cell_cell_ref()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&CELL\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.cell_ref.status == True:
            self.cell_ref.to_input(fout)
        fout.write("\t\t&END CELL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CELL_REF":
                self.cell_ref.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_forcefield_nonbonded_genpot:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&GENPOT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END GENPOT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_forcefield_nonbonded_goodwin:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&GOODWIN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END GOODWIN\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_forcefield_nonbonded_lennard_jones:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LENNARD-JONES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END LENNARD-JONES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_forcefield_nonbonded_williams:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WILLIAMS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END WILLIAMS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_qmmm_forcefield_nonbonded:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.genpot = cp2k_qmmm_forcefield_nonbonded_genpot()
        self.goodwin = cp2k_qmmm_forcefield_nonbonded_goodwin()
        self.lennard_jones = cp2k_qmmm_forcefield_nonbonded_lennard_jones()
        self.williams = cp2k_qmmm_forcefield_nonbonded_williams()

        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&NONBONDED\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.genpot.status == True:
            self.genpot.to_input(fout)
        if self.goodwin.status == True:
            self.goodwin.to_input(fout)
        if self.lennard_jones.status == True:
            self.lennard_jones.to_input(fout)
        if self.williams.status == True:
            self.williams.to_input(fout)
        fout.write("\t\t\t&END NONBONDED\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "GENPOT":
                self.genpot.set_params({item: parmas[item]})
            elif item.split("-")[3] == "GOODWIN":
                self.goodwin.set_params({item: parmas[item]})
            elif item.split("-")[3] == "LENNARD_JONES":
                self.lennard_jones.set_params({item: parmas[item]})
            elif item.split("-")[3] == "WILLIAMS":
                self.williams.set_params({item: parmas[item]})
            else:
                pass




class cp2k_qmmm_forcefield_nonbonded14_genpot:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&GENPOT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END GENPOT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_forcefield_nonbonded14_goodwin:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&GOODWIN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END GOODWIN\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_forcefield_nonbonded14_lennard_jones:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LENNARD-JONES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END LENNARD-JONES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_forcefield_nonbonded14_williams:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WILLIAMS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END WILLIAMS\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_qmmm_forcefield_nonbonded14:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.genpot = cp2k_qmmm_forcefield_nonbonded14_genpot()
        self.goodwin = cp2k_qmmm_forcefield_nonbonded14_goodwin()
        self.lennard_jones = cp2k_qmmm_forcefield_nonbonded14_lennard_jones()
        self.williams = cp2k_qmmm_forcefield_nonbonded14_williams()

        # basic setting


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&NONBONDED14\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.genpot.status == True:
            self.genpot.to_input(fout)
        if self.goodwin.status == True:
            self.goodwin.to_input(fout)
        if self.lennard_jones.status == True:
            self.lennard_jones.to_input(fout)
        if self.williams.status == True:
            self.williams.to_input(fout)
        fout.write("\t\t\t&END NONBONDED14\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "GENPOT":
                self.genpot.set_params({item: parmas[item]})
            elif item.split("-")[3] == "GOODWIN":
                self.goodwin.set_params({item: parmas[item]})
            elif item.split("-")[3] == "LENNARD_JONES":
                self.lennard_jones.set_params({item: parmas[item]})
            elif item.split("-")[3] == "WILLIAMS":
                self.williams.set_params({item: parmas[item]})
            else:
                pass



class cp2k_qmmm_forcefield:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.nonbonded = cp2k_qmmm_forcefiled_nonbonded()
        self.nonbonded14 = cp2k_qmmm_forcefield_nonbonded14()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&FORCEFIELD\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.nonbonded.status == True:
            self.nonbonded.to_input(fout)
        if self.nonbonded14.status == True:
            self.nonbonded14.to_input(fout)
        fout.write("\t\t&END FORCEFIELD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "NONBONDED":
                self.nonbonded.set_params({item: params[item]})
            elif item.split("-")[2] == "NONBONDED14":
                self.nonbonded14.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_force_mixing_buffer_links_link_add_mm_charge:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&ADD_MM_CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END ADD_MM_CHARGE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_force_mixing_buffer_links_link_move_mm_charge:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&MOVE_MM_CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END MOVE_MM_CHARGE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_force_mixing_buffer_links_link:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.add_mm_charge = cp2k_qmmm_force_mixing_buffer_links_link_add_mm_charge()
        self.move_mm_charge = cp2k_qmmm_force_mixing_buffer_links_link_move_mm_charge()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LINK\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.add_mm_charge.status == True:
            self.add_mm_charge.to_input(fout)
        if self.move_mm_charge.status == True:
            self.move_mm_charge.to_input(fout)
        fout.write("\t\t\t\t&END LINK\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "ADD_MM_CHARGE":
                self.add_mm_charge.set_params({item: params[item]})
            elif item.split("-")[4] == "MOVE_MM_CHARGE":
                self.move_mm_charge.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_force_mixing_buffer_links:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.link = cp2k_qmmm_force_mixing_buffer_links_link()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&BUFFER_LINKS\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.link.status == True:
            self.link.to_input(fout)
        fout.write("\t\t\t&END BUFFER_LINKS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_force_mixing_buffer_non_adaptive_link_add_mm_charge:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&ADD_MM_CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END ADD_MM_CHARGE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_force_mixing_buffer_non_adaptive_link_move_mm_charge:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&MOVE_MM_CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END MOVE_MM_CHARGE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_force_mixing_buffer_non_adaptive_link:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.add_mm_charge = cp2k_qmmm_force_mixing_buffer_non_adaptive_link_add_mm_charge()
        self.move_mm_charge = cp2k_qmmm_force_mixing_buffer_non_adaptive_link_move_mm_charge()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&LINK\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.add_mm_charge.status == True:
            self.add_mm_charge.to_input(fout)
        if self.move_mm_charge.status == True:
            self.move_mm_charge.to_input(fout)
        fout.write("\t\t\t\t&END LINK\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "ADD_MM_CHARGE":
                self.add_mm_charge.set_params({item: params[item]})
            elif item.split("-")[4] == "MOVE_MM_CHARGE":
                self.move_mm_charge.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_force_mixing_buffer_non_adaptive_qm_kind:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&QM_KIND\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END QM_KIND\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_force_mixing_buffer_non_adaptive:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.link = cp2k_qmmm_force_mixing_buffer_non_adaptive_link()
        self.qm_kind = cp2k_qmmm_force_mixing_buffer_non_adaptive_qm_kind()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&BUFFER_NON_ADAPTIVE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.link.status == True:
            self.link.to_input(fout)
        if self.qm_kind.status == True:
            self.qm_kind.to_input(fout)
        fout.write("\t\t\t&END BUFFER_NON_ADAPTIVE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "LINK":
                self.link.set_params({item: params[item]})
            elif item.split("-")[3] == "QM_KIND":
                self.qm_kind.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_force_mixing_print_neighbor_lists_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_force_mixing_print_neighbor_lists:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_force_mixing_print_neighbor_lists_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&NEIGHBOR_LISTS\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END NEIGHBOR_LISTS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_force_mixing_print_subcell_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_force_mixing_print_subcell:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_force_mixing_print_subcell_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&SUBCELL\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END SUBCELL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_force_mixing_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.neighbor_lists = cp2k_qmmm_force_mixing_print_neighbor_lists()
        self.subcell = cp2k_qmmm_forcer_mixing_print_subcell()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.neighbor_lists.status == True:
            self.neighbor_lists.to_input(fout)
        if self.subcell.status == True:
            self.subcell.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "NEIGHBOR_LISTS":
                self.neighbor_lists.set_params({item: params[item]})
            elif item.split("-")[3] == "SUBCELL":
                self.subcell.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_force_mixing_qm_non_adaptive_qm_kind:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&QM_KIND\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END QM_KIND\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_force_mixing_qm_non_adaptive:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.qm_kind = cp2k_qmmm_force_mixing_qm_non_adaptive_qm_kind()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&QM_NON_ADAPTIVE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.qm_kind.status == True:
            self.qm_kind.to_input(fout)
        fout.write("\t\t\t&END QM_NON_ADAPTIVE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "QM_KIND":
                self.qm_kind.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_force_mixing_restart_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&RESTART_INFO\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END RESTART_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_qmmm_force_mixing:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.buffer_links = cp2k_qmmm_force_mixing_buffer_links()
        self.buffer_non_adaptive = cp2k_qmmm_force_mixing_buffer_non_adaptive()
        self.printout = cp2k_qmmm_force_mixing_print()
        self.qm_non_adaptive = cp2k_qmmm_force_mixing_qm_non_adaptive()
        self.restart_info = cp2k_qmmm_force_mixing_restart_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&FORCE_MIXING\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.buffer_links.status == True:
            self.buffer_links.to_input(fout)
        if self.buffer_non_adaptive.status == True:
            self.buffer_non_adaptive.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.qm_non_adaptive.status == True:
            self.qm_non_adaptive.to_input(fout)
        if self.restart_info.status == True:
            self.restart_info.to_input(fout)
        fout.write("\t\t&END FORCE_MIXING\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "BUFFER_LINKS":
                self.buffer_links.set_params({item: params[item]})
            elif item.split("-")[2] == "BUFFER_NON_ADAPTIVE":
                self.buffer_non_adaptive.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[2] == "QM_NON_ADAPTIVE":
                self.qm_non_adaptive.set_params({item: params[item]})
            elif item.split("-")[2] == "RESTART_INFO":
                self.restart_info.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_image_charge_eri_mme_cutoff_calib:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CUTOFF_CALIB\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END CUTOFF_CALIB\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_image_charge_eri_mme_eri_mme_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_image_charge_eri_mme_eri_mme_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_image_charge_eri_mme_eri_mme_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CUTOFF_CALIB\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END CUTOFF_CALIB\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_image_charge_eri_mme:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cutoff_calib = cp2k_qmmm_image_charge_eri_mme_cutoff_calib()
        self.eri_mme_info = cp2k_qmmm_image_charge_eri_mme_eri_mme_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ERI_MME\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.cutoff_calib.status == True:
            self.cutoff_calib.to_input(fout)
        if self.eri_mme_info.status == True:
            self.eri_mme_info.to_input(fout)
        fout.write("\t\t\t&END ERI_MME\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CUTOFF_CALIB":
                self.cutoff_calib.set_params({item: params[item]})
            elif item.split("-")[3] == "ERI_MME_INFO":
                self.eri_mme_info.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_image_charge:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.eri_mme = cp2k_qmmm_image_charge_eri_mme()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&IMAGE_CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.eri_mme.status == True:
            self.eri_mme.to_input(fout)
        fout.write("\t\t&END IMAGE_CHARGE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "ERI_MME":
                self.eri_mme.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_interpolator_conv_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_qmmm_interpolator_conv_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_interpolator_conv_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_interpolator_spl_coeffs_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_qmmm_interpolator_spl_coeffs:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_interpolator_spl_coeffs_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SPL_COEFFS\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END SPL_COEFFS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_interpolator:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_qmmm_interpolator_conv_info()
        self.spl_coeffs = cp2k_qmmm_interpolator_spl_coeffs()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        if self.spl_coeffs.status == True:
            self.spl_coeffs.to_input(fout)
        fout.write("\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            elif item.split("-")[2] == "SPL_COEFFS":
                self.spl_coeffs.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_link_add_mm_charge:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&ADD_MM_CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END ADD_MM_CHARGE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_link_move_mm_charge:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MOVE_MM_CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END MOVE_MM_CHARGES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_link:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.add_mm_charge = cp2k_qmmm_link_add_mm_charge()
        self.move_mm_charge = cp2k_qmmm_link_move_mm_charge()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&LINK\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.add_mm_charge.status == True:
            self.add_mm_charge.to_input(fout)
        if self.mvoe_mm_charge.status == True:
            self.move_mm_charge.to_input(fout)
        fout.write("\t\t&END LINK\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "ADD_MM_CHARGE":
                self.add_mm_charge.set_params({item: params[item]})
            elif item.split("-")[2] == "MVOE_MM_CHARGE":
                self.move_mm_charge.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_mm_kind:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MM_KIND\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END MM_KIND\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_check_spline_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_periodic_check_spline:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_periodic_check_spline_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CHECK_SPLINE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END CHECK_SPLINE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_periodic_interpolator_conv_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_interpolator_conv_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_periodic_interpolator_conv_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_periodic_interpolator:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_qmmm_periodic_interpolator_conv_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_periodic_multipole_check_spline_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_periodic_multipole_check_spline:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_periodic_multipole_check_spline_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CHECK_SPLINE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END CHECK_SPLINES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass






class cp2k_qmmm_periodic_multipole_interpolator_conv_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
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


class cp2k_qmmm_periodic_multipole_interpolator_conv_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_periodic_multipole_interpolator_conv_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_qmmm_periodic_multipole_interpolator:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_qmmm_periodic_multipole_interpolator_conv_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.conv_info.status == True:
            self.covn_info.to_input(fout)
        fout.write("\t\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CONV_INFO":
                self.covn_info.set_params({item: params[item]})
            else:
                pass



class cp2k_qmmm_periodic_multipole_program_run_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_periodic_multipole_program_run_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_qmmm_periodic_multipole_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&program_run_info\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&end program_run_info\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_periodic_multipole:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.check_spline = cp2k_qmmm_periodic_multipole_check_splien()
        self.interpolator = cp2k_qmmm_periodic_multipole_interpolator()
        self.program_run_info = cp2k_qmmm_periodic_multipole_program_run_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MULTIPOLE\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.check_spline.status == True:
            self.check_spline.to_input(fout)
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t&END MULTIPOLE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CHECK_SPLINE":
                self.check_spline.set_params({item: params[item]})
            elif item.split("-")[3] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_periodic_poisson_ewald_multipoles:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&MULTIPOLES\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END MULTIPOLES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_poisson_ewald_print_program_run_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_poisson_ewald_print_program_run_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.each = cp2k_qmmm_periodic_poisson_ewald_print_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_periodic_poisson_ewald_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.program_run_info = cp2k_qmmm_periodic_poisson_ewald_print_program_run_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_periodic_poisson_ewald_rs_grid:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&RS_GRID\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t\t&END RS_GRID\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_qmmm_periodic_poisson_ewald:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.multipoles = cp2k_qmmm_periodic_poisson_ewald_multipoles()
        self.printout = cp2k_qmmm_periodic_poisson_ewald_print()
        self.rs_grid = cp2k_qmmm_periodic_poisson_ewald_rs_grid()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EWALD\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.multipoles.status == True:
            self.multipole.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.rs_grid.status == True:
            self.rs_grid.to_input(fout)
        fout.write("\t\t\t\t&END_EALD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "MULTIPOLES":
                self.multipoles.set_params({item: params[item]})
            elif item.split("-")[4] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[4] == "RS_GRID":
                self.rs_grid.set_params({item: params[item]})
            else:
                pass



class cp2k_qmmm_periodic_poisson_implicit_dielectric_dielec_aa_cuboidal:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&DIELEC_AA_CUBOIDAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END DIELEC_AA_CUBOIDAL\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_poisson_implicit_dielectric_dielec_xaa_annular:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&DIELEC_XAA_ANNULAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END DIELEC_XAA_ANNULAR\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_periodic_poisson_implicit_dielectric:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.dielec_aa_cuboidal = cp2k_qmmm_periodic_poisson_implicit_dielectric_dielec_aa__cuboidal()
        self.dielec_xaa_annular = cp2k_qmmm_periodic_poisson_implicit_dielectric_dielec_xaa_annular()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&DIELECTRIC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.dielec_aa_cuboidal.status == True:
            self.dielec_aa_cuboidal.to_input(fout)
        if self.dielec_xaa_annular.status == True:
            self.dielec_xaa_annular.to_input(fout)
        fout.write("\t\t\t\t\t&END DIELECTRIC\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "DIELEC_AA_CUBOIDAL":
                self.dielec_aa_cuboidal.set_params({item: params[item]})
            elif item.split("-")[5] == "DIELEC_XAA_ANNULAR":
                self.dielec_xaa_annular.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_periodic_poisson_implicit_dielectric_bc_aa_cuboidal:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&AA_CUBOIDAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END AA_CUBOIDAL\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_poisson_implicit_dielectric_bc_aa_cylindrical:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&AA_CYLINDRICAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END AA_CYLINDRICAL\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_poisson_implicit_dielectric_bc_aa_planar:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&AA_PLANAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END AA_PLANAR\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_poisson_implicit_dielectric_bc_planar:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&PLANAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END PLANAR\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_periodic_poisson_implicit_dielectric_bc:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.aa_cuboidal = cp2k_qmmm_periodic_poisson_implicit_dielectric_bc_aa_cuboidal()
        self.aa_cylindrical = cp2k_qmmm_periodic_poisson_implicit_dielectric_bc_aa_cylindrical()
        self.aa_planar = cp2k_qmmm_periodic_poisson_implicit_dielectric_bc_aa_planar()
        self.planar = cp2k_qmmm_periodic_poisson_implicit_dielectric_bc_planar()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&DIELECTRIC_BC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.aa_cuboidal.status == True:
            self.aa_cuboidal.to_input(fout)
        if self.aa_cylindrical.status == True:
            self.aa_cylindrical.to_input(fout)
        if self.aa_planar.status == True:
            self.aa_planar.to_input(fout)
        if self.planar.status == True:
            self.planar.to_input(fout)
        fout.write("\t\t\t\t\t&END DIELECTRIC_BC\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "AA_CUBOIDAL":
                self.aa_cuboidal.set_params({item: params[item]})
            elif item.split("-")[5] == "AA_CYLINDRICAL":
                self.aa_cylindrical.set_params({item: params[item]})
            elif item.split("-")[5] == "AA_PLANAR":
                self.aa_planar.set_params({item: params[item]})
            elif item.split("-")[5] == "PLANAR":
                self.planar.set_params({item: params[item]})
            else:
                pass



class cp2k_qmmm_periodic_poisson_implicit:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false

        self.dielectric = cp2k_qmmm_periodic_poisson_implicit_dielectric()
        self.dielectric_bc == cp2k_qmmm_periodic_poisson_implicit_dielectric_bc()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&IMPLICIT\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.dielectric.status == True:
            self.dielectric.to_input(fout)
        if self.dielectric_bc.status == True:
            self.dielectric_bc.to_input(fout)
        fout.write("\t\t\t\t&END IMPLICIT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "DIELECTRIC":
                self.dielectric.set_params({item: params[item]})
            elif item.split("-")[4] == "DIELECTRIC_BC":
                self.dielectric_bc.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_periodic_poisson_mt:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MT\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END MT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass












class cp2k_qmmm_periodic_poisson_multipole_check_spline_each:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_poisson_multipole_check_spline:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.each =cp2k_qmmm_periodic_poisson_multipole_check_spline_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&CHECK_SPLINE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CHECK_SPLINE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_periodic_poisson_multipole_interpolator_conv_info_each:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 8:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_poisson_multipole_interpolator_conv_info:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.each = cp2k_qmmm_periodic_poisson_multipole_interpolator_conv_info_each()
        
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[6] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_periodic_poisson_multipole_interpolator:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.conv_info = cp2k_qmmm_periodic_poisson_multipole_interpolator_conv_info()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_periodic_poisson_multipole_program_run_info_each:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm_periodic_poisson_multipole_program_run_info:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.each = cp2k_qmmm_periodic_poisson_multipole_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_periodic_poisson_multipole:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.check_spline = cp2k_qmmm_periodic_poisson_multipole_check_spline()
        self.interpolator = cp2k_qmmm_periodic_poisson_multipole_interpolator()
        self.program_run_info = cp2k_qmmm_periodic_poisson_multipole_program_run_info()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MULTIPOLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.check_spline.status == True:
            self.check_spline.to_input(fout)
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t\t&END MULTIPOLE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CHECK_SPLINE":
                self.check_spline.set_params({item: params[item]})
            elif item.split("-")[4] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[4] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass



class cp2k_qmmm_periodic_poisson_wavelet:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = false


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&WAVELET\n")
        for item in self.params:
            if self.params[item] is not none:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END WAVELET\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_periodic_poisson:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.ewald = cp2k_qmmm_periodic_poisson_ewald()
        self.implicit = cp2k_qmmm_periodic_poisson_implicit()
        self.mt = cp2k_qmmm_periodic_poisson_mt()
        self.multipole = cp2k_qmmm_periodic_poisson_multipole()
        self.wavelet = cp2k_qmmm_periodic_poisson_wavelet()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&POISSON\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.ewald.status == True:
            self.ewald.to_input(fout)
        if self.implicit.status == True:
            self.implicit.to_input(fout)
        if self.mt.status == True:
            self.mt.to_input(fout)
        if self.multipole.status == True:
            self.multipole.to_input(fout)
        if self.wavelet.status == True:
            self.wavelet.to_input(fout)
        fout.write("\t\t\t&END POISSON\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EWALD":
                self.ewald.set_params({item: params[item]})
            elif item.split("-")[3] == "IMPLICIT":
                self.implicit.set_params({item: params[item]})
            elif item.split("-")[3] == "MT":
                self.mt.set_params({item: params[item]})
            elif item.split("-")[3] == "MULTIPOLE":
                slef.multipole.set_params({item: params[item]})
            elif item.split("-")[3] == "WAVELET":
                self.wavelet.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_periodic:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.check_spline = cp2k_qmmm_periodic_check_spline()
        self.interpolator = cp2k_qmmm_periodic_interpolator()
        self.multipole = cp2k_qmmm_periodic_multipole()
        self.poisson = cp2k_qmmm_periodic_poisson()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PERIODIC\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.check_spline.status == True:
            self.check_spline.to_input(fout)
        if self.interpolator.status == True:
            selfinterpolator.to_input(fout)
        if self.multipole.status == True:
            self.multipole.to_input(fout)
        if self.poisson.status == True:
            self.poisson.to_input(fout)
        fout.write("\t\t&END PERIODIC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CHECK_SPLINE":
                self.check_spline.set_params({item: params[item]})
            elif item.split("-")[2] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[2] == "MULTIPOLE":
                self.multipole.set_params({item: params[item]})
            elif item.split("-")[2] == "POISSON":
                self.poisson.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_derivatives_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_derivatives:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_derivatives_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_print_dipole_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_dipole:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_dipole_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_print_grid_information_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_grid_information:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_grid_information_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_image_charge_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_image_charge_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_image_charge_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_image_charge_restart_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_image_charge_restart:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_image_charge_restart_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_mm_potential_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_mm_potential:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_mm_potential_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_periodic_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_periodic_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_periodic_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_pgf_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_pgf:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_pgf_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_potential_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_potential:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_potential_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_program_banner_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_program_banner:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_program_banner_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_program_run_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_program_run_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_qmmm_charges_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_qmmm_charges:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_qmmm_charges_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_qmmm_link_info_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_qmmm_link_info:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_qmmm_link_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_qmmm_matrix_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_qmmm_matrix:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_qmmm_matrix_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_print_qs_derivatives_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_print_qs_derivatives:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_qmmm_print_qs_derivatives_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_qmmm_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.derivatives = cp2k_qmmm_print_derivatives()
        self.dipole = cp2k_qmmm_print_dipole()
        self.grid_information = cp2k_qmmm_print_grid_information()
        self.image_charge_info = cp2k_qmmm_print_image_charge_info()
        self.image_charge_restart = cp2k_qmmm_print_image_charge_restart()
        self.mm_potential = cp2k_qmmm_print_mm_potential()
        self.periodic_info = cp2k_qmmm_print_periodic_info()
        self.pgf = cp2k_qmmm_print_pgf()
        self.potential = cp2k_qmmm_print_potential()
        self.program_banner = cp2k_qmmm_print_program_banner()
        self.program_run_info = cp2k_qmmm_print_program_run_info()
        self.qmmm_charges = cp2k_qmmm_print_qmmm_charges()
        self.qmmm_link_info = cp2k_qmmm_print_qmmm_link_info()
        self.qmmm_matrix = cp2k_qmmm_print_qmmm_matrix()
        self.qs_derivatives = cp2k_qmmm_print_qs_derivatives()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.derivatives.status == True:
            self.drivatives.to_input(fout)
        if self.dipole.status == True:
            self.dipole.to_input(fout)
        if self.grid_information.status == True:
            self.grid_information.to_input(fout)
        if self.image_charge_info.status == True:
            self.image_charge_info.to_input(fout)
        if self.image_charge_restart.status == True:
            self.image_charge_restart.to_input(fout)
        if self.mm_potential.status == True:
            self.mm_potential.to_input(fout)
        if self.periodic_info.status == True:
            self.periodic_info.to_input(fout)
        if self.pgf.status == True:
            self.pgg.to_input(fout)
        if self.potential.status == True:
            self.potential.to_input(fout)
        if self.program_banner.status == True:
            self.program_banner.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.qmmm_charges.status == True:
            self.qmmm_charges.to_input(fout)
        if self.qmmm_link_info.status == True:
            self.qmmm_link_info.to_input(fout)
        if self.qmmm_matrix.status == True:
            self.qmmm_matrix.to_input(fout)
        if self.qs_derivatives.status == True:
            self.qs_derivatives.to_input(fout)
        fout.write("\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "DERIVATIVES":
                self.drivatives.set_params({item: params[item]})
            elif item.split("-")[2] == "DIPOLE":
                self.dipole.set_params({item: params[item]})
            elif item.split("-")[2] == "GRID_INFORMATION":
                self.grid_information.set_params({item: params[item]})
            elif item.split("-")[2] == "IMAGE_CHARGE_INFO":
                self.image_charge_info.set_params({item: params[item]})
            elif item.split("-")[2] == "IMAGE_CHARGE_RESTART":
                self.image_charge_restart.set_params({item: params[item]})
            elif item.split("-")[2] == "MM_POTENTIAL":
                self.mm_potential.set_params({item: params[item]})
            elif item.split("-")[2] == "PERIODIC_INFO":
                self.periodic_info.set_params({item: params[item]})
            elif item.split("-")[2] == "PGF":
                self.pgf.set_params({item: params[item]})
            elif item.split("-")[2] == "POTENTIAL":
                self.potential.set_params({item: params[item]})
            elif item.split("-")[2] == "PROGRAM_BANNER":
                self.program_banner.set_params({item: params[item]})
            elif item.split("-")[2] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[2] == "QMMM_CHARGES":
                self.qmmm_charges.set_params({item: params[item]})
            elif item.split("-")[2] == "QMMM_LINK_INFO":
                self.qmmm_link_info.set_params({item: params[item]})
            elif item.split("-")[2] == "QMMM_MATRIX":
                self.qmmm_matrix.set_params({item: params[item]})
            elif item.split("-")[2] == "QS_DERIVATIVES":
                self.qs_derivatives.set_params({item: params[item]})
            else:
                pass

class cp2k_qmmm_qm_kind:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&QM_KIND\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END QM_KIND\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_qmmm_walls:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False


        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&WALLS\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END WALLS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_qmmm:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.cell = cp2k_qmmm_cell()
        self.forcefield = cp2k_qmmm_forcefield()
        self.force_mixing = cp2k_qmmm_force_mixing()
        self.image_charge = cp2k_qmmm_image_charge()
        self.interpolator = cp2k_qmmm_interpolator()
        self.link = cp2k_qmmm_link()
        self.mm_kind = cp2k_qmmm_mm_kind()
        self.periodic = cp2k_qmmm_periodic()
        self.printout = cp2k_qmmm_print()
        self.qm_kind = cp2k_qmmm_qm_kind()
        self.walls = cp2k_qmmm_walls()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&QMMM\n")
        for item in self.params:
            if self.params[item] is not None:                    
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.cell.status == True:
            self.cell.to_input(fout)
        if self.forcefield.status == True:
            self.forcefield.to_input(fout)
        if self.force_mixing.status == True:
            self.force_mixing.to_input(fout)
        if self.image_charge.status == True:
            self.image_charge.to_input(fout)
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.link.status == True:
            self.link.to_input(fout)
        if self.mm_kind.status == True:
            self.mm_kind.to_input(fout)
        if self.periodic.status == True:
            self.periodic.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.qm_kind.status == True:
            self.qm_kind.to_input(fout)
        if self.walls.status == True:
            self.walls.to_input(fout)
        fout.write("\t&END QMMM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "CELL":
                self.cell.set_params({item: params[item]})
            elif item.split("-")[1] == "FORCEFIELD":
                self.forcefield.set_params({item: params[item]})
            elif item.split("-")[1] == "FORCE_MIXING":
                self.force_mixing.set_params({item: params[item]})
            elif item.split("-")[1] == "IMAGE_CHARGE":
                self.image_charge.set_params({item: params[item]})
            elif item.split("-")[1] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[1] == "LINK":
                self.link.set_params({item: params[item]})
            elif item.split("-")[1] == "MM_KIND":
                self.mm_kind.set_params({item: params[item]})
            elif item.split("-")[1] == "PERIODIC":
                self.periodic.set_params({item: params[item]})
            elif item.split("-")[1] == "QM_KIND":
                self.qm_kind.set_params({item: params[item]})
            elif item.split("-")[1] == "WALLS":
                self.walls.set_params({item: params[item]})
            else:
                pass
