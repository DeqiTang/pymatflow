#!/usr/bin/env python
# _*_ coding: utf-8 _*_

# =================================
# CP2K / FORCE_EVAL / DFT / POISSON
# =================================


class cp2k_dft_poisson_ewald_multipoles:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&MULTIPOLES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END MULTIPOLES\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_poisson_ewald_print_program_run_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_poisson_ewald_print_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_poisson_ewald_print_program_run_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass


class cp2k_dft_poisson_ewald_print:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.program_run_info = cp2k_dft_poisson_ewald_print_program_run_info()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t\t&END PRINT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_poisson_ewald_rs_grid:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&RS_GRID\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t&END RS_GRID\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_poisson_ewald:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.multipoles = cp2k_dft_poisson_ewald_multipoles()
        self.printout = cp2k_dft_poisson_ewald_print()
        self.rs_grid = cp2k_dft_poisson_ewald_rs_grid()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&EWALD\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.multipoles.status == True:
            self.multipole.to_input(fout)
        if self.printout.status == True:
            self.printout.to_input(fout)
        if self.rs_grid.status == True:
            self.rs_grid.to_input(fout)
        fout.write("\t\t\t&END EWALD\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "MULTIPOLES":
                self.multipoles.set_params({item: params[item]})
            elif item.split("-")[3] == "PRINT":
                self.printout.set_params({item: params[item]})
            elif item.split("-")[3] == "RS_GRID":
                self.rs_grid.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_poisson_implicit_dielectric_dielec_aa_cuboidal:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&DIELEC_AA_CUBOIDAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END DIELEC_AA_CUBOIDAL\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_poisson_implicit_dielectric_dielec_xaa_annular:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&DIELEC_XAA_ANNULAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END DIELEC_XAA_ANNULAR\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_poisson_implicit_dielectric:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.dielec_aa_cuboidal = cp2k_dft_poisson_implicit_dielectric_dielec_aa_cuboidal()
        self.dielec_xaa_annular = cp2k_dft_poisson_implicit_dielectric_dielec_xaa_annular()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DIELECTRIC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.dielec_aa_cuboidal.status == True:
            self.dielec_aa_cuboidal.to_input(fout)
        if self.dielec_xaa_annular.status == True:
            self.dielec_xaa_annular.to_input(fout)
        fout.write("\t\t\t\t&END DIELECTRIC\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "DIELEC_AA_CUBOIDAL":
                self.dielec_aa_cuboidal.set_params({item: params[item]})
            elif item.split("-")[4] == "DIELEC_XAA_ANNULAR":
                self.dielec_xaa_annular.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_poisson_implicit_dirichlet_bc_aa_cuboidal:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&AA_CUBOIDAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END AA_CUBOIDAL\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_poisson_implicit_dirichlet_bc_aa_cylindrical:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&AA_CYLINDRICAL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END AA_CYLINDRICAL\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_poisson_implicit_dirichlet_bc_aa_planar:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&AA_PLANAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END AA_PLANAR\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_poisson_implicit_dirichlet_bc_planar:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&PLANAR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t&END PLANAR\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_poisson_implicit_dirichlet_bc:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.aa_cuboidal = cp2k_dft_poisson_implicit_dirichlet_bc_aa_cuboidal()
        self.aa_cylindrical = cp2k_dft_poisson_implicit_dirichlet_bc_aa_cylindrical()
        self.aa_planar = cp2k_dft_poisson_implicit_dirichlet_bc_aa_planar()
        self.planar = cp2k_dft_poisson_implicit_dirichlet_bc_planar()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&DIRICHLET_BC\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.aa_cuboidal.status == True:
            self.aa_cuboidal.to_input(fout)
        if self.aa_cylindrical.status == True:
            self.aa_cylindrical.to_input(fou)
        if self.aa_planar.status == True:
            self.aa_planar.to_input(fout)
        if self.planar.status == True:
            self.planar.to_input(fout)
        fout.write("\t\t\t\t&END DIRICHLET_BC\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "AA_CUBOIDAL":
                self.cuboidal.set_params({item: params[item]})
            elif item.split("-")[4] == "AA_CYLINDRICAL":
                self.aa_cylindrical.set_params({item: params[item]})
            elif item.split("-")[4] == "AA_PLANAR":
                self.aa_planar.set_params({item: params[item]})
            elif item.split("-")[4] == "PLANAR":
                self.planar.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_poisson_implicit:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.dielectric = cp2k_dft_poisson_implicit_dielectric()
        self.dirichlet_bc = cp2k_dft_poisson_implicit_dirichlet_bc()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&IMPLICIT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.dielectric.status == True:
            self.dielectric.to_input(fout)
        if self.dirichlet_bc.status == True:
            self.dirichlet_bc.to_input(fout)
        fout.write("\t\t\t&END IMPLICIT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "DIELECTRIC":
                self.dielectric.set_params({item: params[item]})
            elif item.split("-")[3] == "DIRICHLET_BC":
                self.dirichlet_bc.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_poisson_mt:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END MT\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_poisson_multipole_check_spline_each:
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

class cp2k_dft_poisson_multipole_check_spline:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_poisson_multipole_check_spline_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&CHECK_SPLINE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END CHECK_SPLINE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_poisson_multipole_interpolator_conv_info_each:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t\t&EACH\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t\t\t\t&END EACH\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 7:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_poisson_multipole_interpolator_conv_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_poisson_multipole_interpolator_conv_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t\t&CONV_INVO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t\t&END CONV_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 6:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[5] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_poisson_multipole_interpolator:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.conv_info = cp2k_dft_poisson_multipole_interpolator_conv_info()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&INTERPOLATOR\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.conv_info.status == True:
            self.conv_info.to_input(fout)
        fout.write("\t\t\t\t&END INTERPOLATOR\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "CONV_INFO":
                self.conv_info.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_poisson_multipole_program_run_info_each:
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

class cp2k_dft_poisson_multipole_program_run_info:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.each = cp2k_dft_poisson_multipole_program_run_info_each()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t\t&PROGRAM_RUN_INFO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.each.status == True:
            self.each.to_input(fout)
        fout.write("\t\t\t\t&END PROGRAM_RUN_INFO\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 5:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[4] == "EACH":
                self.each.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_poisson_multipole:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.chekc_spline = cp2k_dft_poisson_multipole_check_spline()
        self.interpolator = cp2k_dft_poisson_multipole_interpolator()
        self.program_run_info = cp2k_dft_poisson_multipole_program_run_info()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&MULTIPOLE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.check_spline.status == True:
            self.check_spline.to_input(fou)
        if self.interpolator.status == True:
            self.interpolator.to_input(fout)
        if self.program_run_info.status == True:
            self.program_run_info.to_input(fout)
        fout.write("\t\t\t&END MULTIPOLE\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[3] == "CHECK_SPLINE":
                self.check_spline.set_params({item: params[item]})
            elif item.split("-")[3] == "INTERPOLATOR":
                self.interpolator.set_params({item: params[item]})
            elif item.split("-")[3] == "PROGRAM_RUN_INFO":
                self.program_run_info.set_params({item: params[item]})
            else:
                pass

class cp2k_dft_poisson_wavelet:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t\t&WAVELET\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, str(self.params[item])))
        fout.write("\t\t\t&END WAVELET\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_poisson:
    def __init__(self):
        self.params = {
                "PERIODIC": "XYZ",
                "POISSON_SOLVER": "PERIODIC",
                }
        self.status = False

        self.ewald = cp2k_dft_poisson_ewald()
        self.implicit = cp2k_dft_poisson_implicit()
        self.mt = cp2k_dft_poisson_mt()
        self.multipole = cp2k_dft_poisson_multipole()
        self.wavelet = cp2k_dft_poisson_wavelet()

        # basic setting

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&POISSON\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, str(self.params[item])))
        if self.ewald.status == True:
            self.ewald.to_input(fout)
        if self.implicit.status == True:
            self.implicit.to_input(fout)
        if self.mt.status == True:
            self.mt.to_input(fout)
        if self.multipole.status == True:
            self.multipole.to_input(fout)
        if self.wavelet.status == True:
            self.wavelet.to_input(fout)
        fout.write("\t\t&END POISSON\n")

    def set_params(self, params):
        #
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "EWALD":
                self.ewald.set_params({item: params[item]})
            elif item.split("-")[2] == "IMPLICIT":
                self.implicit.set_params({item: params[item]})
            elif item.split("-")[2] == "MT":
                self.mt.set_params({item: params[item]})
            elif item.split("-")[2] == "MULTIPOLE":
                self.multipole.set_params({item: params[item]})
            elif item.split("-")[2] == "WAVELET":
                self.wavelet.set_params({item: params[item]})
            else:
                pass
