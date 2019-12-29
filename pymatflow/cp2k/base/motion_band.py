#!/usr/bin/env python
# _*_ coding: utf-8 _*_

from pymatflow.base.xyz import base_xyz

# ====================
# CP2K / MOTION / BAND
# ====================

class cp2k_motion_band_banner_each:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not none:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_band_banner:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_band_banner_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&BANNER\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END BANNER\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_band_ci_neb:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&CI_NEB\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END CI_NEB\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_band_convergence_control:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&CONVERGENCE_CONTROL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END CONVERGENCE_CONTROL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_band_convergence_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_band_convergence_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_band_convergence_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&CONVERGENCE_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END CONVERGENCE_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_band_energy_each:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_band_energy:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_band_energy_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&ENERGY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END ENERGY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_band_optimize_band_diis_diis_info_each:
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
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_band_optimize_band_diis_diis_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_band_optimize_band_diis_diis_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DIIS_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END DIIS_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_band_optimize_band_diis:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.diis_info = cp2k_motion_band_optimize_band_diis_diis_info()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&DIIS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.diis_info.status == True:
            self.diis_info.to_input(fout)
        fout.write("\t\t\t&END DIIS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "DIIS_INFO":
                self.diis_info.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_band_optimize_band_md_temp_control:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&TEMP_CONTROL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END TEMP_CONTROL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_band_optimize_band_md_vel_control:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&VEL_CONTROL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END VEL_CONTROL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_band_optimize_band_md:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.temp_control = cp2k_motion_band_optimize_band_md_temp_control()
        self.vel_control = cp2k_motion_band_optimize_band_md_vel_control()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.temp_control.status == True:
            self.temp_contro.to_input(fout)
        if self.vel_control.status == True:
            self.vel_control.to_input(fout)
        fout.write("\t\t\t&END MD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "TEMP_CONTROL":
                self.temp_control.set_params({item: params[item]})
            elif item.split("-")[2] == "VEL_CONTROL":
                self.vel_control.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_band_optimize_band:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.diis = cp2k_motion_band_optimize_band_diis()
        self.md = cp2k_motion_band_optimize_band_md()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&OPTIMIZE_BAND\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.diis.status == True:
            self.diis.to_input(fout)
        if self.md.status == True:
            self.md.to_input(fout)
        fout.write("\t\t&END OPTIMIZE_BAND\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "DIIS":
                self.diis.set_params({item: params[item]})
            elif item.split("-")[1] == "MD":
                self.md.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_band_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_band_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_band_program_run_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_motion_band_replica_coord:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&COORD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END COORD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_band_replica_velocity:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&VELOCITY\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END VELOCITY\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_motion_band_replica:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.coord = cp2k_motion_band_replica_coord()
        self.velocity = cp2k_motion_band_replica_velocity()
        # basic setting
        
    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&REPLICA\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.coord.status == True:
            self.coord.to_input(fout)
        if self.velocity.status == True:
            self.velocity.to_input(fout)
        fout.write("\t\t&END REPLICA\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "COORD":
                self.coord.set_params({item: params[item]})
            elif item.split("-")[1] == "VELOCITY":
                self.velocity.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_band_replica_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END EACH\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_band_replica_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_motion_band_replica_info_each()
        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&REPLICA_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t&END REPLICA_INFO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_motion_band_string_method:
    def __init__(self):
        self.params = {
                }
        self.status = False


    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&STRING_METHOD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t&END STRING_METHOD\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_motion_band:
    def __init__(self):
        self.params = {
                "ALIGN_FRAMES": None,
                "BAND_TYPE": None, # CI-NEB, IT-NEB, SM
                "K_SPRING": None,
                "NPROC_REP": None,
                "NUMBER_OF_REPLICA": None,
                "POT_TYPE": None,
                "PROC_DIST_TYPE": None,
                "ROTATE_FRAMES": None,
                "USE_COLVARS": None,
                }
        self.status = False
        
        self.banner = cp2k_motion_band_banner()
        self.ci_neb = cp2k_motion_band_ci_neb()
        self.convergence_control = cp2k_motion_band_convergence_control()
        self.convergence_info = cp2k_motion_band_convergence_info()
        self.energy = cp2k_motion_band_energy()
        self.optimize_band = cp2k_motion_band_optimize_band()
        self.program_run_info = cp2k_motion_band_program_run_info()
        self.replica = cp2k_motion_band_replica()
        self.replica_info = cp2k_motion_band_replica_info()
        self.string_method = cp2k_motion_band_string_method()

        self.images = []

        self.params["BAND_TYPE"] = "CI-NEB"
        self.params["NUMBER_OF_REPLICA"] = 10
        self.params["K_SPRING"] = 0.05
        self.params["ROTATE_FRAMES"] = "TRUE"
        self.params["ALIGN_FRAMES"] = "TRUE"

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t&BAND\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t%s %s\n" % (item, str(self.params[item])))
        if self.banner.status == True:
            self.banner.to_input(fout)
        if self.ci_neb.status == True:
            self.ci_neb.to_input(fout)
        if self.convergence_control.status == True:
            self.convergence_control.to_input(fout)
        if self.convergence_info.status == True:
            self.convergence_info.to_input(fout)
        if self.energy.status == True:
            self.energy.to_input(fout)
        if self.optimize_band.status == True:
            self.optimize_band.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        for image in self.images:
            fout.write("\t\t&REPLICA\n")
            fout.write("\t\t\tCOORD_FILE_NAME %s\n" % image.file)
            fout.write("\t\t&END REPLICA\n")
        if self.replica_info.status == True:
            self.replica_info.to_input(fout)
        if self.string_method.status == True:
            self.string_method.to_input(fout)
        fout.write("\t&END BAND\n")

    def get_images(self, images):
        for image in images:
            xyz = base_xyz()
            xyz.get_xyz(image)
            self.images.append(xyz)

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 2:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[1] == "BANNER":
                self.banner.set_params({item: params[item]})
            elif item.split("-")[1] == "CI_NEB":
                self.ci_neb.set_params({item: params[item]})
            elif item.split("-")[1] == "CONVERGENCE_CONTROL":
                self.convergence_control.set_params({item: params[item]})
            elif item.split("-")[1] == "CONVERGENCE_INFO":
                self.convergence_info.set_params({item: params[item]})
            elif item.split("-")[1] == "ENERGY":
                self.energy.set_params({item: params[item]})
            elif item.split("-")[1] == "OPTIMIZE_BAND":
                self.optimize_band.set_params({item: params[item]})
            elif item.split("-")[1] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            elif item.split("-")[1] == "REPLICA":
                self.replica.set_params({item: params[item]})
            elif item.split("-")[1] == "REPLICA_INFO":
                self.replica_info.set_params({item: params[item]})
            elif item.split("-")[1] == "STRING_METHOD":
                self.string_method.set_params({item: params[item]})
            else:
                pass


