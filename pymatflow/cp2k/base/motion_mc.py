#!/usr/bin/env python
# _*_ coding: utf-8 _*_



class cp2k_motion_mc_avbmc:
    def __init__(self):
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&AVBMC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END AVBMC\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_mc_max_displacements_box_probabilities:
    def __init__(self):
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&BOX_PROBABILITIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END BOX_PROBABILITIES\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_mc_max_displacements_mol_probabilities:
    def __init__(self):
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MOL_PROBABILITIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END MOL_PROBABILITIES\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_mc_max_displacements:
    def __init__(self):
        self.status = False

        self.box_probabilities = cp2k_motion_mc_max_displacements_box_probabilities()
        self.mol_Probabilities = cp2k_motion_mc_max_displacements_mol_probabilities()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MAX_DISPLACEMENTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.box_probabilities.status == True:
            self.box_probabilities.to_input(fout)
        if self.mol_probabilities.status == True:
            self.mol_probabilities.to_input(fout)
        fout.write("\t\t&END MAX_DISPLACEMENTS\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "BOX_PROBABILITIES":
                self.box_probabilities.set_params({item: params[item]})
            elif item.split("-")[2] == "MOL_PROBABILITIES":
                self.mol_probabilities.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_mc_move_probabilities_box_probabilities:
    def __init__(self):
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&BOX_PROBABILITIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END BOX_PROBABILITIES\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_mc_move_probabilities_mol_probabilities:
    def __init__(self):
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MOL_PROBABILITIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END MOL_PROBABILITIES\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_mc_move_probabilities:
    def __init__(self):
        self.status = False

        self.box_probabilities = cp2k_motion_mc_move_probabilities_box_probabilities()
        self.mol_Probabilities = cp2k_motion_mc_move_probabilities_mol_probabilities()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MOVE_PROBABILITIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.box_probabilities.status == True:
            self.box_probabilities.to_input(fout)
        if self.mol_probabilities.status == True:
            self.mol_probabilities.to_input(fout)
        fout.write("\t\t&END MOVE_PROBABILITIES\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "BOX_PROBABILITIES":
                self.box_probabilities.set_params({item: params[item]})
            elif item.split("-")[2] == "MOL_PROBABILITIES":
                self.mol_probabilities.set_params({item: params[item]})
            else:
                pass



class cp2k_motion_mc_move_updates:
    def __init__(self):
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&MOVE_UPDATES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END MOVE_UPDATES\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_mc:
    def __init__(self):
        self.status = False

        self.avmbc = cp2k_motion_mc_avbmc()
        self.max_displacements = cp2k_motion_mc_max_displacements()
        self.move_probabilities = cp2k_motion_mc_move_probabilities()
        self.move_updates = cp2k_motion_mc_move_updates()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&MC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        if self.avbmc.status == True:
            self.avbmc.to_input(fout)
        if self.max_displacements.status == True:
            self.max_displacements.to_input(fout)
        if self.move_probabilities.status == True:
            self.move_probabilities.to_input(fout)
        if self.move_updates.status == True:
            self.move_update.to_input(fout)
        fout.write("\t&END MC\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "AVBMC":
                self.avbmc.set_params({item: params[item]})
            elif item.split("-")[1] == "MAX_DISPLACEMENTS":
                self.max_displacements.set_params({item: params[item]})
            elif item.split("-")[1] == "MOVE_PROBABILITIES":
                self.move_probabilities.set_params({item: params[item]})
            elif item.split("-")[1] == "MOVE_UPDATE":
                self.move_updates.set_params({item: params[item]})
            else:
                pass

