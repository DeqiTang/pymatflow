#!/usr/bin/env python
# _*_ coding: utf-8 _*_

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

    def to_properties(self, fout):
        fout.write("\t\t&RESP\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&SLAB_SAMPLING\n")
        fout.write("\t\t\t\tRANGE 0.3 3.0\n")
        fout.write("\t\t\t\tATOM_LIST 2 3 5 8\n")
        fout.write("\t\t\t\tSURF_DIRECTION Z\n")
        fout.write("\t\t\t&END SLAB_SAMPLING\n")

        fout.write("\t\t\t&SLAB_SAMPLING\n")
        fout.write("\t\t\t\tRANGE 1.0 3.0\n")
        fout.write("\t\t\t\tATOM_LIST 2 3 5 8\n")
        fout.write("\t\t\t\tSURF_DIRECTION -Z\n")
        fout.write("\t\t\t&END SLAB_SAMPLING\n")

        fout.write("\t\t\t&PRINT\n")
        fout.write("\t\t\t\t&V_RESP_CUBE\n")
        fout.write("\t\t\t\t&END V_RESP_CUBE\n")
        fout.write("\t\t\t&END PRINT\n")
        fout.write("\t\t&END RESP\n")


class cp2k_properties_linres:
    def __init__(self):
        self.params = {
                "ENERGY_GAP": None,
                "EPS": None,
                "MAX_ITER": None,
                "PRECONDITIONER": None,
                }
        self.status = False

    def to_properties(self, fout):
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
        self.resp = cp2k_properties_resp()
        self.linres = cp2k_properties_linres()

    def to_force_eval(self, fout):
        fout.write("\t&PROPERTIES\n")
        if self.resp.status == True:
            self.resp.to_properties(fout)
        if self.linres.status == True:
            self.linres.to_properties(fout)
        fout.write("\t&END PROPERTIES\n")
