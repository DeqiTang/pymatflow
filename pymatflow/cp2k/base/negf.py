#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil

"""
Usage:
"""

class cp2k_negf_contact_bulk_region:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&BULK_REGION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END BULK_REGION\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_negf_contact_print_dos_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_negf_contact_print_dos:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_negf_contact_print_dos_each()
        # basic setting
    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&DOS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status ==True:
            self.each.to_input(fout)
        fout.write("\t\t\t&END DOS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_negf_contact_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.dos = cp2k_negf_contact_print_dos()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.dos.status == True:
            self.dos.to_input(fout)
        fout.write("\t\t&END PRINT\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "DOS":
                self.dos.set_params({item: params[item]})
            else:
                pass

class cp2k_negf_contact_screening_region:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&SCREENING_REGION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END SCREENING_REGION\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_negf_contact:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.bulk_region = cp2k_negf_contact_bulk_region()
        self.printout = cp2k_negf_contact_print()
        self.screening_region = cp2k_negf_contact_screening_region()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&CONTACT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.bulk_region.status == True:
            self.bulk_region.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.screening_region.status == True:
            self.screening_region.to_input(fout)
        fout.write("\t&END CONTACT\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "BULK_REGION":
                self.bulk_region.set_params({item: params[item]})
            elif item.split("-")[1] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[1] == "SCREENING_REGION":
                self.screening_region.set_params({item: params[item]})
            else:
                pass


class cp2k_negf_mixing:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&MIXING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END MIXING\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_negf_print_dos_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_negf_print_dos:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_negf_print_dos_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&DOS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END DOS\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_negf_print_transmission_each:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EACH\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_negf_print_transmission:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_negf_print_transmission_each()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t\t&TRANSMISSION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END TRANSMISSION\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_negf_print:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.dos = cp2k_negf_print_dos()
        self.transmission = cp2k_negf_print_transmission()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.dos.status == True:
            self.dos.to_input(fout)
        if self.transmission.status == True:
            self.transmission.to_input(fout)
        fout.write("\t&END PRINT\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "DOS":
                self.dos.set_params({item: params[item]})
            elif item.split("-")[1] == "TRANSMISSION":
                self.transmission.set_params({item: params[item]})
            else:
                pass


class cp2k_negf_scattering_region:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("\t&SCATTERING_REGION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t&END SCATTERING_REGION\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_negf:
    """

    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.contact = cp2k_negf_contact()
        self.mixing = cp2k_negf_mixing()
        self.printout = cp2k_negf_print()
        self.scattering_region = cp2k_negf_scattering_region()
        # basic setting

    def to_input(self, fout):
        # fout: a file stream for writing
        fout.write("&NEGF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t%s %s\n" % (item, self.params[item]))
        if self.contact.status == True:
            self.contact.to_input(fout)
        if self.mixing.status == True:
            self.mixing.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.scattering_region.status == True:
            self.scattering_region.to_input(fout)
        fout.write("&END NEGF\n")
        fout.write("\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 1:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[0] == "CONTACT":
                self.contact.set_params({item: params[item]})
            elif item.split("-")[0] == "MIXING":
                self.mixing.set_params({item: params[item]})
            elif ietm.split("-")[0] == "PRITN":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[0] == "SCATTERING_REGION":
                self.scattering_region.set_params({item: params[item]})
            else:
                pass


