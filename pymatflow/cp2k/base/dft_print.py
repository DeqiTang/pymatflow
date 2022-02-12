#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np

# ===============================
# CP2K / FORCE_EVAL / DFT / PRINT
# ===============================

class cp2k_dft_print_active_space:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&ACTIVE_SPACE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&ACTIVE_SPACE\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_band_structure:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.params["ADDED_MOS"] = 10
        self.params["FILE_NAME"] = "bands.bs"

    def to_input(self, fout):
        # so when you activate band calculation by setting self.status to True,
        # you should also pass an kpath to self.set_band()
        # to set the kpoints path used in band structure calculation 
        fout.write("\t\t\t&BAND_STRUCTURE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        for i in range(len(self.kpath) - 1):
            if self.kpath[i][4] != "|" and type(self.kpath[i][4] == int):
                fout.write("\t\t\t\t&KPOINT_SET\n")
                fout.write("\t\t\t\t\tUNITS B_VECTOR\n")
                fout.write("\t\t\t\t\tSPECIAL_POINT %s %f %f %f\n" % (self.kpath[i][3], self.kpath[i][0], self.kpath[i][1], self.kpath[i][2]))
                fout.write("\t\t\t\t\tSPECIAL_POINT %s %f %f %f\n" % (self.kpath[i+1][3], self.kpath[i+1][0], self.kpath[i+1][1], self.kpath[i+1][2]))
                fout.write("\t\t\t\t\tNPOINTS %d\n" % (self.kpath[i][4]-1))
                fout.write("\t\t\t\t&END KPOINT_SET\n")
        fout.write("\t\t\t&END BAND_STRUCTURE\n")

    def set_band(self, kpath=None):
        """
        kpath: 
            the high symmetry k point path used in bands structure calculation
            in format like this:
            
            [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]
            
            if connect_indicator in a kpoint is an integer, then it will connect to the following point
            through the number of kpoints defined by connect_indicator.
            
            if connect_indicator in a kpoint is '|', then it will not connect to the following point,
        """
        self.band_structure = True
        self.kpath = kpath

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_efield_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&EFIELD_CUBE ON\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EFIELD_CUBE\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_elf_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

        # basic setting
        self.params["STRIDE"] = [1, 1, 1]

    def to_input(self, fout):
        fout.write("\t\t\t&ELF_CUBE ON\n")
        for item in self.params:
            if self.params[item] is not None:
                if item == "STRIDE":
                    fout.write("\t\t\t\t%s %d %d %d\n" % (item, self.params[item][0], self.params[item][1], self.params[item][2]))
                else:
                    fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END ELF_CUBE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_e_density_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.params["STRIDE"] = [1, 1, 1]
        self.params["FILENAME"] = None
        self.params["APPEND"] = "TRUE"

    def to_input(self, fout):
        fout.write("\t\t\t&E_DENSITY_CUBE ON\n")
        for item in self.params:
            if self.params[item] is not None:
                if item == "STRIDE":
                    fout.write("\t\t\t\t%s %d %d %d\n" % (item, self.params[item][0], self.params[item][1], self.params[item][2]))
                else:
                    fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END E_DENSITY_CUBE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_pdos:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.params["NLUMO"] = -1 # print all the projected DOS available
        self.params["COMPONENTS"] = "" # split the density by quantum number

    def to_input(self, fout):
        fout.write("\t\t\t&PDOS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END PDOS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_moments:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.params["PERIODIC"] = "FALSE"

    def to_input(self, fout):
        fout.write("\t\t\t&MOMENTS\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END MOMENTS\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_electric_field_gradient:
    def __init__(self):
        self.params = {
                }
        self.status = False
        
        self.params["FILENAME"] = "electric-field-gradient.dat"

    def to_input(self, fout):
        fout.write("\t\t\t&ELECTRIC_FIELD_GRADIENT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END ELECTRIC_FIELD_GRADIENT\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_print_external_potential_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.params["FILENAME"] ="external-potential.cube"
        
    def to_input(self, fout):
        fout.write("\t\t\t&EXTERNAL_POTENTIAL_CUBE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END EXTERNAL_POTENTIAL_CUBE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_print_lowdin:
    def __init__(self):
        self.params = {
                }
        self.status = False

        
    def to_input(self, fout):
        fout.write("\t\t\t&LOWDIN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END LOWDIN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_mo:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.params["FILENAME"] ="mo.dat"
        
    def to_input(self, fout):
        fout.write("\t\t\t&MO\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END MO\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass

class cp2k_dft_print_mo_cubes:
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.params["NHOMO"] = -1
        self.params["NLUMO"] = -1
        
    def to_input(self, fout):
        fout.write("\t\t\t&MO_CUBES\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END MO_CUBES\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_mulliken:
    def __init__(self):
        self.params = {
                }
        self.status = False

        
    def to_input(self, fout):
        fout.write("\t\t\t&MULLIKEN\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END MULLIKEN\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_stm:
    def __init__(self):
        self.params = {
                }
        self.status = False

        
    def to_input(self, fout):
        fout.write("\t\t\t&STM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END STM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_subcell:
    def __init__(self):
        self.params = {
                }
        self.status = False

        
    def to_input(self, fout):
        fout.write("\t\t\t&SUBCELL\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END SUBCELL\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_dft_print_tot_density_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

        
    def to_input(self, fout):
        fout.write("\t\t\t&TOT_DENSITY_CUBE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END TOT_DENSITY_CUBE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_dft_print_v_hartree_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

        
    def to_input(self, fout):
        fout.write("\t\t\t&V_HARTREE_CUBE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END V_HARTREE_CUBE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_v_xc_cube:
    def __init__(self):
        self.params = {
                }
        self.status = False

        
    def to_input(self, fout):
        fout.write("\t\t\t&V_XC_CUBE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END V_XC_CUBE\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_wannier90:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&WANNIER90\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END WANNIER90\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass


class cp2k_dft_print_wfn_mix:
    def __init__(self):
        self.params = {
                }
        self.status = False

    def to_input(self, fout):
        fout.write("\t\t\t&WFN_MIX\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END WFN_MIX\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_dft_print_xray_diffraction_spectrum:
    def __init__(self):
        self.params = {
                }
        self.status = False

        
    def to_input(self, fout):
        fout.write("\t\t\t&XRAY_DIFFRACTION_SPECTRUM\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t&END XRAY_DIFFRACTION_SPECTRUM\n")

    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 4:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass



class cp2k_dft_print:
    """
    analysis toos:
        cubecruncher:
            the cubecruncher is a tool to do various operations on cube files.
            eg. it can help with analysis of charge density cube files.
            it lies in cp2k/tool/cubecruncher
    """
    def __init__(self):
        self.params = {
                }
        self.status = False

        self.active_space = cp2k_dft_print_active_space()

        self.adjmat_write = False
        self.band_structure = cp2k_dft_print_band_structure()
        self.basis_molopt_quantities = False
        self.basis_set_file = False
        self.derivatives = False
        self.dft_control_parameters = False

        self.efield_cube = cp2k_dft_print_efield_cube()

        self.electric_field_gradient = cp2k_dft_print_electric_field_gradient()

        self.elf_cube = cp2k_dft_print_elf_cube()

        self.energy_windows = False

        self.external_potential_cube = cp2k_dft_print_external_potential_cube()

        self.e_density_cube = cp2k_dft_print_e_density_cube()

        self.gapw = False
        self.hirshfeld = False
        self.hyperfine_coupling_tensor = False
        self.implicit_psolver = False
        self.kinetic_energy = False
        self.kpoints = False
        self.ks_csr_write = False

        self.lowdin = cp2k_dft_print_lowdin()

        self.mao_analysis = False
        self.minbas_analysis = False

        self.mo = cp2k_dft_print_mo()

        self.moments = cp2k_dft_print_moments()

        self.mo_cubes = cp2k_dft_print_mo_cubes()

        self.mulliken = cp2k_dft_print_mulliken()

        self.neighbor_lists = False
        self.optimize_lri_basis = False
        self.overlap_condition = False
        
        self.pdos = cp2k_dft_print_pdos()
        
        self.plus_u = False
        self.program_banner = False
        self.sccs = False

        self.stm = cp2k_dft_print_stm()

        self.subcell = cp2k_dft_print_subcell()

        self.s_csr_write = False

        self.tot_density_cube = cp2k_dft_print_tot_density_cube()

        self.v_hartree_cube = cp2k_dft_print_v_hartree_cube()

        self.v_xc_cube = cp2k_dft_print_v_xc_cube()

        self.wannier90 = cp2k_dft_print_wannier90()

        self.wfn_mix = cp2k_dft_print_wfn_mix()

        self.xray_diffraction_spectrum = cp2k_dft_print_xray_diffraction_spectrum()

    def to_input(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))

        if self.active_space.status == True:
            self.active_space.to_input(fout)

        if self.band_structure.status == True:
            self.band_structure.to_input(fout)
        
        if self.efield_cube.status == True:
            self.efield_cube.to_input(fout)

        if self.electric_field_gradient.status == True:
            self.electric_field_gradient.to_input(fout)

        if self.elf_cube.status == True:
            self.elf_cube.to_input(fout)

        if self.e_density_cube.status == True:
            self.e_density_cube.to_input(fout)

        if self.external_potential_cube.status == True:
            self.external_potential_cube.to_input(fout)

        if self.pdos.status == True:
            self.pdos.to_input(fout)

        if self.lowdin.status == True:
            self.lowdin.to_input(fout)

        if self.mo.status == True:
            self.mo.to_input(fout)

        if self.moments.status == True:
            self.moments.to_input(fout)

        if self.mo_cubes.status == True:
            self.mo_cubes.to_input(fout)

        if self.mulliken.status == True:
            self.mulliken.to_input(fout)

        if self.stm.status == True:
            self.stm.to_input(fout)

        if self.subcell.status == True:
            self.subcell.to_input(fout)
       
        if self.tot_density_cube.status == True:
            self.tot_density_cube.to_input(fout)

        if self.v_hartree_cube.status == True:
            self.v_hartree_cube.to_input(fout)

        if self.v_xc_cube.status == True:
            self.v_xc_cube.to_input(fout)

        if self.wannier90.status == True:
            self.wannier90.to_input(fout)

        if self.wfn_mix.status == True:
            self.wfn_mix.to_input(fout)

        if self.xray_diffraction_spectrum.status == True:
            self.xray_diffraction_spectrum.to_input(fout)

        fout.write("\t\t&END PRINT\n")


    def set_params(self, params):
        for item in params:
            if len(item.split("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "BAND_STRUCTURE":
                self.band_structure.set_params({item: params[item]})
            elif item.split("-")[2] == "ACTIVE_SPACE":
                self.active_space.set_params({item: params[item]})
            elif item.split("-")[2] == "EFIELD_CUBE":
                self.efield_cube.set_params({item: params[item]})
            elif item.split("-")[2] == "ELF_CUBE":
                self.elf_cube.set_params({item: params[item]})
            elif item.split("-")[2] == "E_DENSITY_CUBE":
                self.e_density_cube.set_params({item: params[item]})
            elif item.split("-")[2] == "LOWDIN":
                self.lowdin.set_params({item: params[item]})
            elif item.split("-")[2] == "MO":
                self.mo.set_params({item: params[item]})           
            elif item.split("-")[2] == "MO_CUBES":
                self.mo_cubes.set_params({item: params[item]})           
            elif item.split("-")[2] == "MOMENTS":
                self.moments.set_params({item: params[item]})           
            elif item.split("-")[2] == "MULLIKEN":
                self.mulliken.set_params({item: params[item]})           
            elif item.split("-")[2] == "ELECTRIC_FIELD_GRADIENT":
                self.electric_field_gradient.set_params({item: params[item]})           
            elif item.split("-")[2] == "EXTERNAL_POTENTIAL_CUBE":
                self.external_potential_cube.set_params({item: params[item]})           
            elif item.split("-")[2] == "PDOS":
                self.pdos.set_params({item: params[item]})           
            elif item.split("-")[2] == "SUBCELL":
                self.subcell.set_params({item: params[item]})           
            elif item.split("-")[2] == "STM":
                self.stm.set_params({item: params[item]})           
            elif item.split("-")[2] == "TOT_DENSITY-CUBE":
                self.tot_density_cube.set_params({item: params[item]})           
            elif item.split("-")[2] == "V_HARTREE_CUBE":
                self.v_hartree_cube.set_params({item: params[item]})           
            elif item.split("-")[2] == "V_XC_CUBE":
                self.x_xc_cube.set_params({item: params[item]})           
            elif item.split("-")[2] == "WANNIER90":
                self.wannier90.set_params({item: params[item]})           
            elif item.split("-")[2] == "WFN_MIX":
                self.wfn_mix.set_params({item: params[item]})           
            elif item.split("-")[2] == "XRAY_DIFFRACTION_SPECTRUM":
                self.xray_diffraction_spectrum.set_params({item: params[item]})           
            else:
                pass

