#!/usr/bin/env python
# _*_ coding: utf-8 _*_



class cp2k_motion_md_thermostat:
    def __init__(self):
        self.params = {
                "REGION": None,
                "TYPE": None,
                }
        self.params["TYPE"] = "NOSE"

    def to_motion_md(self, fout):
        fout.write("\t\t&THERMOSTAT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t&END THERMOSTAT\n")


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
                "STEPS": 50,
                "STEP_START_VAL": None,
                "TEMPERATURE": None,
                "TEMPERATURE_ANNEALING": None,
                "TEMP_KIND": None,
                "TEMP_TOL": None,
                "TIMESTEP": None,
                "TIME_START_VAL": None,
                }
        self.thermostat = cp2k_motion_md_thermostat()

    def to_motion(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&MD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        if self.params["ENSEMBLE"] == "NVT":
            self.thermostat.to_motion_md(fout)
        fout.write("\t&END MD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]

