#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# ===================
# ===================

class cp2k_mm_forcefield_bend_ub:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&UB\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END UB\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mm_forcefield_bend:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.ub = cp2k_mm_forcefield_bend_ub()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&BEND\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.ub.status == True:
            self.ub.to_input(fout)
        fout.write("\t\t\t&END BEND\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "UB":
                self.ub.set_params({item: params[item]})
            else:
                pass


class cp2k_mm_forcefield_bond:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&BOND\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END BOND\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_charge:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CHARGE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END CHARGE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_charges:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&CHARGES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END CHARGES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mm_forcefield_dipole_damping:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DAMPING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END DAMPING\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mm_forcefield_dipole:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.damping = cp2k_mm_forcefield_dipole_damping()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&DIPOLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.damping.status == True:
            self.damping.to_input(fout)
        fout.write("\t\t\t&END DIPOLE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DAMPING":
                self.damping.set_params({item: params[item]})
            else:
                pass

class cp2k_mm_forcefield_improper:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&IMPROPER\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END IMPROPER\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mm_forcefield_nonbonded_bmhft:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&BMHFT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END BMHFT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mm_forcefield_nonbonded_bmhftd:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&BMHFTD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END BMHFTD\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mm_forcefield_nonbonded_buck4ranges:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&BUCK4RANGES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END BUCK4RANGES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mm_forcefield_nonbonded_buckmorse:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&BUCKMORSE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END BUCKMORSE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_nonbonded_eam:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&EAM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END EAM\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_nonbonded_genpot:
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

class cp2k_mm_forcefield_nonbonded_goodwin:
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

class cp2k_mm_forcefield_nonbonded_ipbv:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&IPBV\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END IPBV\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_mm_forcefield_nonbonded_lennard_jones:
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

class cp2k_mm_forcefield_nonbonded_quip:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&QUIP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END QUIP\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_nonbonded_siepmann:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&SIEPMANN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END SIEPMANN\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_nonbonded_tersoff:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&TERSOFF\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END TERSOFF\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_nonbonded_williams:
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

class cp2k_mm_forcefield_nonbonded:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.bmhft = cp2k_mm_forcefield_nonbonded_bmhft()
        self.bmhftd = cp2k_mm_forcefield_nonbonded_bmhftd()
        self.buck4ranges = cp2k_mm_forcefield_nonbonded_buck4ranges()
        self.buckmorse = cp2k_mm_forcefield_nonbonded_buckmorse()
        self.eam = cp2k_mm_forcefield_nonbonded_eam()
        self.genpot = cp2k_mm_forcefield_nonbonded_genpot()
        self.goodwin = cp2k_mm_forcefield_nonbonded_goodwin()
        self.ipbv = cp2k_mm_forcefield_nonbonded_ipbv()
        self.lennard_jones = cp2k_mm_forcefield_nonbonded_lennard_jones()
        self.quip = cp2k_mm_forcefield_nonbonded_quip()
        self.siepmann = cp2k_mm_forcefield_nonbonded_siepmann()
        self.tersoff = cp2k_mm_forcefield_nonbonded_tersoff()
        self.williams = cp2k_mm_forcefield_nonbonded_williams()

        # basic setting
        

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&NONBONDED\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.bmhft.status == True:
            self.bmhft.to_input(fout)
        if self.bmhftd.status == True:
            self.bmhftd.to_input(fout)
        if self.buck4ranges.status == True:
            self.buck4ranges.to_input(fout)
        if self.buckmorse.status == True:
            self.buckmorse.to_input(fout)
        if self.eam.status == True:
            self.eam.to_input(fout)
        if self.genpot.status == True:
            self.genpot.to_input(fout)
        if self.goodwin.status == True:
            self.goodwin.to_input(fout)
        if self.ipbv.status == True:
            self.ipbv.to_input(fout)
        if self.lennard_jones.status == True:
            self.lennard_jones.to_input(fout)
        if self.quip.status == True:
            self.quip.to_input(fout)
        if self.siepmann.status == True:
            self.siepmann.to_input(fout)
        if self.tersoff.status == True:
            self.tersoff.to_input(fout)
        if self.williams.status == True:
            self.williams.to_input(fout)
        fout.write("\t\t\tEND NONBONDED\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "BMHFT":
                self.bmhft.set_params({item: parmas[item]})
            elif item.split("-")[3] == "BMHFTD":
                self.bmhftd.set_params({item: parmas[item]})
            elif item.split("-")[3] == "BUCK4RANGES":
                self.buck4ranges.set_params({item: parmas[item]})
            elif item.split("-")[3] == "BUCKMORSE":
                self.buckmorse.set_params({item: parmas[item]})
            elif item.split("-")[3] == "EAM":
                self.eam.set_params({item: parmas[item]})
            elif item.split("-")[3] == "GENPOT":
                self.genpot.set_params({item: parmas[item]})
            elif item.split("-")[3] == "GOODWIN":
                self.goodwin.set_params({item: parmas[item]})
            elif item.split("-")[3] == "IPBV":
                self.ipbv.set_params({item: parmas[item]})
            elif item.split("-")[3] == "LENNARD_JONES":
                self.lennard_jones.set_params({item: parmas[item]})
            elif item.split("-")[3] == "QUIP":
                self.quip.set_params({item: parmas[item]})
            elif item.split("-")[3] == "SIEPMANN":
                self.siepmann.set_params({item: parmas[item]})
            elif item.split("-")[3] == "TERSOFF":
                self.tersoff.set_params({item: parmas[item]})
            elif item.split("-")[3] == "WILLIAMS":
                self.williams.set_params({item: parmas[item]})
            else:
                pass



class cp2k_mm_forcefield_nonbonded14_genpot:
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

class cp2k_mm_forcefield_nonbonded14_goodwin:
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


class cp2k_mm_forcefield_nonbonded14_lennard_jones:
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


class cp2k_mm_forcefield_nonbonded14_williams:
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



class cp2k_mm_forcefield_nonbonded14:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.genpot = cp2k_mm_forcefield_nonbonded14_genpot()
        self.goodwin = cp2k_mm_forcefield_nonbonded14_goodwin()
        self.lennard_jones = cp2k_mm_forcefield_nonbonded14_lennard_jones()
        self.williams = cp2k_mm_forcefield_nonbonded14_williams()

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

                




class cp2k_mm_forcefield_opbend:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&OPBEND\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END OPBEND\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_quadrupole:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&QUADRUPOLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END QUADRUPOLE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_shell:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SHELL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SHELL\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_spline:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&SPLINE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SPLINE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield_torsion:
    def __init__(self):
        self.params = {
                }       
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&TORSION\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END TORSION\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_mm_forcefield:
    def __init__(self):
        self.params = {
                }       
        self.status = False

        self.bend = cp2k_mm_forcefield_bend()
        self.bond = cp2k_mm_forcefield_bond()
        self.charge = cp2k_mm_forcefield_charge()
        self.charges = cp2k_mm_forcefield_charges()
        self.dipole = cp2k_mm_forcefield_dipole()
        self.improper = cp2k_mm_forcefield_improper()
        self.nonbonded = cp2k_mm_forcefield_nonbonded()
        self.nonbonded14 = cp2k_mm_forcefield_nonbonded14()
        self.opbend = cp2k_mm_forcefield_opbend()
        self.quadrupole = cp2k_mm_forcefield_quadrupole()
        self.shell = cp2k_mm_forcefield_shell()
        self.spline = cp2k_mm_forcefield_spline()
        self.torsion = cp2k_mm_forcefield_torsion()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&FORCEFIELD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.bend.status == True:
            self.bend.to_input(fout)
        if self.bond.status == True:
            self.bond.to_input(fout)
        if self.charge.status == True:
            self.charge.to_input(fout)
        if self.charges.status == True:
            self.charges.to_input(fout)
        if self.dipole.status == True:
            self.dipole.to_input(fout)
        if self.improper.status == True:
            self.improper.to_input(fout)
        if self.nonbonded.status == True:
            self.nonbonded.to_input(fout)
        if self.nonbonded14.status == True:
            self.nonbonded14.to_input(fout)
        if self.opbend.status == True:
            self.opbend.to_input(fout)
        if self.quadrupole.status == True:
            self.quadrupole.to_input(fout)
        if self.shell.status == True:
            self.shell.to_input(fout)
        if self.torsion.status == True:
            self.torsion.to_input(fout)
        fout.write("\t\t&END FORCEFIELD\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "BEND":
                self.bend.set_params({item: params[item]})
            elif item.split("-")[2] == "BOND":
                self.bond.set_params({item: params[item]})
            elif item.split("-")[2] == "CHARGE":
                self.charge.set_params({item: params[item]})
            elif item.split("-")[2] == "CHARGES":
                self.charges.set_params({item: params[item]})
            elif item.split("-")[2] == "DIPOLE":
                self.dipole.set_params({item: params[item]})
            elif item.split("-")[2] == "IMPROPER":
                self.improper.set_params({item: params[item]})
            elif item.split("-")[2] == "NONBONDED":
                self.nonbonded.set_params({item: params[item]})
            elif item.split("-")[2] == "NONBONDED14":
                self.nonbonded14.et_params({item: params[item]})
            elif item.split("-")[2] == "OPBEND":
                self.opbend.set_params({item: params[item]})
            elif item.split("-")[2] == "QUADRUPOLE":
                self.quadrupole.set_params({item: params[item]})
            elif item.split("-")[2] == "SHELL":
                self.shell.set_params({item: params[item]})
            elif item.split("-")[2] == "SPLINE":
                self.spline.set_params({item: params[item]})
            elif item.split("-")[2] == "TORSION":
                self.torsion.set_params({item: params[item]})
            else:
                pass
