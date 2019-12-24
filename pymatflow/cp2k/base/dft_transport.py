#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ==============================
# ==============================


class cp2k_dft_transport_beyn:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&BEYN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END BEYN\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_transport_contact:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CONTACT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END CONTACT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_transport_pexsi:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PEXSI\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END PEXSI\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_transport_print_current_each:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_transport_print_current:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.each = cp2k_dft_transport_print_current_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CURRENT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END CURRENT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_transport_print:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.current = cp2k_dft_transport_print_current()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.current.status == True:
            self.current.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CURRENT":
                self.current.set_params({item: params[item]})
            else:
                pass



class cp2k_dft_transport:
    def __init__(self):
        self.params = {
                "TRANSPORT_METHOD": None,
                }       
        self.status = False

        self.beyn = cp2k_dft_transport_beyn()
        self.contact = cp2k_dft_transport_contact()
        self.pexsi = cp2k_dft_transport_pexsi()
        self.printout = cp2k_dft_transport_print()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&TRANSPORT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.beyn.status == True:
            self.beyn.to_input(fout)
        if self.contact.status == True:
            self.contact.to_input(fout)
        if self.pexsi.status == True:
            self.pexsi.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t&END TRANSPORT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "BEYN":
                self.beyn.set_params({item: params[item]})
            elif item.split("-")[2] == "CONTACT":
                self.contact.set_params({item: params[item]})
            elif item.split("-")[2] == "PEXSI":
                self.pexsi.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass
