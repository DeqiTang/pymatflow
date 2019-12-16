#!/usr/bin/env python
# _*_ coding: utf-8 _*_

class cp2k_properties_resp_constraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&CONSTRAINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END CONSTRAINT\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_resp_print_v_resp_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&V_RESP_CUBE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END V_RESP_CUBE\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_resp_print:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.v_resp_cube = cp2k_properties_resp_print_v_resp_cube()

        # basic setting
        self.v_resp_cube.status = True


    def to_input(self, fout):
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.v_resp_cube.status == True:
            self.v_resp_cube.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "V_RESP_CUBE":
                self.v_resp_cube.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_resp_restraint:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&RESTARINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END RESTRAINT\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_resp_slab_sampling:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
        self.params["RANGE"] = [0.3, 3.0]
        self.params["ATOM_LIST"] = [1]  # list like: [1, 2 ,3 , 4]
        self.params["SURF_DIRECTION"] = "Z"


    def to_input(self, fout):
        fout.write("\t\t\t&SLAB_SAMPLING\n")
        for item in self.params:
            if self.params[item] is not None:
                if item == "RANGE":
                    fout.write("\t\t\t\t%s %f %f\n" % (item, self.params[item][0], self.params[item][1]))
                elif item == "ATOM_LIST":
                    fout.write("\t\t\t\t%s " % item)
                    for atom in self.params[item]:
                        fout.write(" %d" % atom)
                    fout.write("\n")
                else:
                    fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SLAB_SAMPLING\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_resp_sphere_sampling:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&SPHERE_SAMPLING\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SPHERE_SAMPLING\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_resp:
    def __init__(self):
        self.params = {
                "INTEGER_TOTAL_CHARGE": None,
                "RESTRAIN_HEAVIES_STRENGTH": None,
                "RESTRAIN_HEAVIES_TO_ZERO": None,
                "STRIDE": None,
                "USE_REPEAT_METHOD": None,
                "WIDTH": None,
                }
        self.status = False

        self.constraint = cp2k_properties_resp_constraint()
        self.printout = cp2k_properties_resp_print()
        self.restraint = cp2k_properties_resp_restraint()
        self.slab_sampling = cp2k_properties_resp_slab_sampling()
        self.sphere_sampling = cp2k_properties_resp_sphere_sampling()

        # basic setting
        self.printout.status = True
        self.slab_sampling.status = True

    def to_input(self, fout):
        fout.write("\t\t&RESP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.constraint.status == True:
            self.constraint.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.restraint.status == True:
            self.restraint.to_input(fout)
        if self.slab_sampling.status == True:
            self.slab_sampling.to_input(fout)
        if self.sphere_sampling.status == True:
            self.sphere_samplint.to_input(fout)
        fout.write("\t\t&END RESP\n")

        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CONSTRAINT":
                self.constraint.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[2] == "RESTRAINT":
                self.restraint.set_params({item: params[item]})
            elif item.split("-")[2] == "SLAB_SAMPLING":
                self.slab_sampling.set_params({item: params[item]})
            elif item.split("-")[2] == "SPHERE_SAMPLING":
                self.sphere_sampling.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres:
    def __init__(self):
        self.params = {
                "ENERGY_GAP": None,
                "EPS": None,
                "MAX_ITER": None,
                "PRECONDITIONER": None,
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t&LINRES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&POLAR\n")
        fout.write("\t\t\t\tDO_RAMAN .TRUE.\n")
        fout.write("\t\t\t&END POLAR\n")

        fout.write("\t\t\t&CURRENT\n")
        fout.write("\t\t\t&END CURRENT\n")

        fout.write("\t\t&END LINRES\n")

class cp2k_properties:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.resp = cp2k_properties_resp()
        self.linres = cp2k_properties_linres()

    def to_input(self, fout):
        fout.write("\t&PROPERTIES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, self.params[item]))
        if self.resp.status == True:
            self.resp.to_input(fout)
        if self.linres.status == True:
            self.linres.to_input(fout)
        fout.write("\t&END PROPERTIES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "RESP":
                self.resp.set_params({item: params[item]})
            elif item.split("-")[1] == "LINRES":
                self.linres.set_params({item: params[item]})
            else:
                pass
