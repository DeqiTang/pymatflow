#!/usr/bin/env python
# _*_ coding: utf-8 _*_


class cp2k_motion_md_thermostat_nose:
    def __init__(self):
        self.params = {
                "LENGTH": None,
                "MTS": None,
                "TIMECON": None,
                "YOSHIDA": None,
                }
        self.status = False
        self.params["LENGTH"] = 3
        self.params["MTS"] = 2
        self.params["TIMECON"] = 1.0e3
        self.params["YOSHIDA"] = 3

    def to_input(self, fout):
        fout.write("\t\t\t&NOSE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END NOSE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "VELOCITY":
                print("=======================================\n")
                print("warning: cp2k_motion_md_thermostat_nose\n")
                print("VELOCITY section not implemented yet\n")
                sys.exit(1)
            elif item.split("-")[3] == "MASS":
                print("=======================================\n")
                print("warning: cp2k_motion_md_thermostat_nose\n")
                print("MASS section not implemented yet\n")
                sys.exit(1)
            elif item.split("-")[3] == "FORCE":
                print("=======================================\n")
                print("warning: cp2k_motion_md_thermostat_nose\n")
                print("FORCE section not implemented yet\n")
                sys.exit(1)
            elif item.split("-")[3] == "COORD":
                print("=======================================\n")
                print("warning: cp2k_motion_md_thermostat_nose\n")
                print("COORD section not implemented yet\n")
                sys.exit(1)

class cp2k_motion_md_thermostat:
    def __init__(self):
        self.params = {
                "REGION": None,
                "TYPE": None,
                }
        self.status = False
        self.params["TYPE"] = "NOSE"
        self.nose = cp2k_motion_md_thermostat_nose()

    def to_input(self, fout):
        fout.write("\t\t&THERMOSTAT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        if self.params["TYPE"].upper() == "NOSE":
            self.nose.to_input(fout)

        fout.write("\t\t&END THERMOSTAT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "NOSE":
                self.nose.set_params({item: params[item]})

class cp2k_motion_md_reftraj:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t&REFTRAJ\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END REFTRAJ\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]

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
        self.status = False
        self.thermostat = cp2k_motion_md_thermostat()
        self.reftraj = cp2k_motion_md_reftraj()
        # basic default setting
        self.params["TIMESTEP"] = 0.5
        self.params["STEPS"] = 1000

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&MD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        if self.params["ENSEMBLE"] == "NVT":
            self.thermostat.to_input(fout)
        if self.params["ENSEMBLE"] == "REFTRAJ":
            self.reftraj.to_input(fout)
        fout.write("\t&END MD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "THERMOSTAT":
                self.thermostat.set_params({item: params[item]})
