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

class cp2k_properties_resp_print_coord_fit_points_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_resp_print_coord_fit_points:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_resp_print_coord_fit_points_each()
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&COORD_FIT_POINTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END COORD_FIT_POINTS\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "COORD_FIT_POINTS":
                self.coord_fit_points.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_resp_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_resp_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_resp_print_program_run_info_each()
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PROGRAM_RUN_INFO\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_resp_print_resp_charges_to_file_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_resp_print_resp_charges_to_file:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_resp_print_resp_charges_to_file_each()
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&RESP_CHARGES_TO_FILE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END RESP_CHARGES_TO_FILE\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_resp_print_v_resp_cube_each:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END EACH\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_resp_print_v_resp_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_resp_print_v_resp_cube_each()
        
        # basic setting


    def to_input(self, fout):
        fout.write("\t\t\t\t&V_RESP_CUBE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END V_RESP_CUBE\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_resp_print:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.coord_fit_points = cp2k_properties_resp_print_coord_fit_points()
        self.program_run_info = cp2k_properties_resp_print_program_run_info()
        self.resp_charges_to_file = cp2k_properties_resp_print_resp_charges_to_file()
        self.v_resp_cube = cp2k_properties_resp_print_v_resp_cube()

        # basic setting
        self.v_resp_cube.status = True


    def to_input(self, fout):
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.coord_fit_points.status == True:
            self.coord_fit_points.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        if self.resp_charges_to_file.status == True:
            self.resp_charges_to_file.to_input(fout)
        if self.v_resp_cube.status == True:
            self.v_resp_cube.to_input(fout)
        fout.write("\t\t\t&END PRINT\n")
        
    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "COORD_FIT_POINTS":
                self.coord_fit_points.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[3] == "RESP_CHARGES_TO_FIEL":
                self.resp_charges_to_file.set_params({item: params[item]})
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

class cp2k_properties_linres_current_interpolator_conv_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_current_interpolator_conv_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_current_interpolator_conv_info_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_current_interpolator:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_properties_linres_current_interpolator_conv_info()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_current_print_current_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_current_print_current_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_current_print_current_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CURRENT_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CURRENT_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_current_print_response_function_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_current_print_response_function_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_current_print_response_function_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&RESPONSE_FUNCTION_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END RESPONSE_FUNCTION_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_current_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.current_cubes = cp2k_properties_linres_current_print_currenct_cubes()
        self.response_function_cubes = cp2k_properties_linres_current_print_response_function_cubes()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.current_cubes.status == True:
            self.current_cubes.to_input(fout)
        if self.response_function_cubes.status == True:
            self.response_function_cubes.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CURRENT_CUBES":
                self.current_cubes.set_params({item: params[item]})
            elif item.split("-")[4] == "RESPONSE_FUNCTION_CUBES":
                self.response_function_cubes.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_current:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.interpolator = cp2k_properties_linres_current_interpolator()
        self.printout = cp2k_properties_linres_current_print()
        
        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&CURRENT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        fout.write("\t\t\t&END CURRENT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_interpolator_conv_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_interpolator_conv_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_interpolator_conv_info_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&CONV_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_epr_interpolator:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_properties_linres_epr_interpolator_conv_info()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_g_tensor_xc:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&XC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END XC\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_g_tensor:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_g_tensor_each()
        self.xc = cp2k_properties_linres_epr_print_g_tensor_xc()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&G_TENSOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        if self.xc.status == True:
            self.xc.to_input(fout)
        fout.write("\t\t\t\t\t&END G_TENSOR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            elif item.split("-")[5] == "XC":
                self.xc.set_params({item: params[item]})
            else:
                pass



class cp2k_properties_linres_epr_print_nablavks_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_properties_linres_epr_print_nablavks_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_nablavks_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&NABLAVKS_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END NABLAVKS_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_properties_linres_epr_print_response_function_cubes_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_epr_print_response_function_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_properties_linres_epr_print_response_function_cubes_each()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t\t&RESPONSE_FUNCTION_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END RESPONSE_FUNCTION_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_properties_linres_epr_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.g_tensor = cp2k_properties_linres_epr_print_g_tensor()
        self.nablavks_cubes = cp2k_properties_linres_epr_print_nablavks_cubes()
        self.response_function_cubes = cp2k_properties_linres_epr_print_response_function_cubes()

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.g_tensor.status == True:
            self.g_tensor.to_input(fout)
        if self.nablavks_cubes.status == True:
            self.nablavks_cubes.to_input(fout)
        if self.response_function_cubes.status == True:
            self.response_function_cubes.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "G_TENSOR":
                self.g_tensor.set_params({item: params[item]})
            elif item.split("-")[4] == "NABLAVKS_CUBES":
                self.nablavks_cubes.set_params({item: params[item]})
            elif item.split("-")[4] == "RESPONSE_FUNCTION_CUBES":
                self.response_function_cubes.set_params({item: params[item]})
            else:
                pass



class cp2k_properties_linres_epr:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&EPR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EPR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_localize:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&LOCALIZE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END LOCALIZE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_nmr:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&NMR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END NMR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_polar:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&POLAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END POLAR\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END PRINT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres_spinspin:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting

    def to_input(self, fout):
        fout.write("\t\t\t&SPINSPIN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END SPINSPIN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_properties_linres:
    def __init__(self):
        self.params = {
                "ENERGY_GAP": none,
                "EPS": none,
                "MAX_ITER": none,
                "PRECONDITIONER": None,
                }
        self.status = False

        self.current = cp2k_properties_linres_current()
        self.epr = cp2k_properties_linres_epr()
        self.loclalize = cp2k_properties_linres_localize()
        self.nmr = cp2k_properties_linres_nmr()
        self.polar = cp2k_properties_linres_polar()
        self.printout = cp2k_properties_linres_print()
        self.spinspin = cp2k_properties_linres_spinspin()

        # basic setting
        self.polar.status = True
        self.current.status = True

    def to_input(self, fout):
        fout.write("\t\t&LINRES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.current.status == True:
            self.current.to_input(fout)
        if self.epr.status == True:
            self.epr.to_input(fout)
        if self.localize.status == True:
            self.localize.to_input(fout)
        if self.nmr.status == True:
            self.nmr.to_input(fout)
        if self.polar.status == True:
            self.polar.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.spinspin.status == True:
            self.spinspin.to_input(fout)
        fout.write("\t\t&END LINRES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "CURRENT":
                self.current.set_params({item: params[item]})
            elif item.split("-")[2] == "EPR":
                self.epr.set_params({item: params[item]})
            elif item.split("-")[2] == "LOCALIZE":
                self.localize.set_params({item: params[item]})
            elif item.split("-")[2] == "NMR":
                self.nmr.set_params({item: params[item]})
            elif item.split("-")[2] == "POLAR":
                self.polar.set_params({item: params[item]})
            elif item.split("-")[2] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[2] == "SPINSPIN":
                self.spinspin.set_params({item: params[item]})
            else:
                pass

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
