#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import numpy as np
import seekpath

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
        # self.kpoints_seekpath is defined by self.set_band()
        # so when you activate band calculation by setting self.status to True,
        # you should also pass an instance of base_xyz to self.set_band()
        # to set the kpoints path using seekpath program
        fout.write("\t\t\t&BAND_STRUCTURE\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t\t%s %s\n" % (item, self.params[item]))
        fout.write("\t\t\t\t&KPOINT_SET\n")
        fout.write("\t\t\t\t\tUNITS B_VECTOR\n")
        point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][0][0]]
        fout.write("\t\t\t\t\tSPECIAL_POINT %s %f %f %f\n" % (self.kpoints_seekpath["path"][0][0], point[0], point[1], point[2]))
        point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][0][1]]
        fout.write("\t\t\t\t\tSPECIAL_POINT %s %f %f %f\n" % (self.kpoints_seekpath["path"][0][1], point[0], point[1], point[2]))
        for i in range(1, len(self.kpoints_seekpath["path"])):
            if self.kpoints_seekpath["path"][i][0] == self.kpoints_seekpath["path"][i-1][1]:
                point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][1]]
                fout.write("\t\t\t\t\tSPECIAL_POINT %s %f %f %f\n" % (self.kpoints_seekpath["path"][i][1], point[0], point[1], point[2]))
            else:
                point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][0]]
                fout.write("\t\t\t\t\tSPECIAL_POINT %s %f %f %f\n" % (self.kpoints_seekpath["path"][i][0], point[0], point[1], point[2]))
                point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][1]]
                fout.write("\t\t\t\t\tSPECIAL_POINT %s %f %f %f\n" % (self.kpoints_seekpath["path"][i][1], point[0], point[1], point[2]))
        fout.write("\t\t\t\t\tNPOINTS 20\n")
        fout.write("\t\t\t\t&END KPOINT_SET\n")
        fout.write("\t\t\t&END BAND_STRUCTURE\n")

    def set_band(self, xyz):
        """
        xyz is an instance of base_xyz, it is used by seekpath to generate k-path
        """
        self.band_structure = True
        # --------------
        # using seekpath
        # --------------
        lattice = [xyz.cell[0:3], xyz.cell[3:6], xyz.cell[6:9]]
        positions = []
        numbers = []
        a = np.sqrt(xyz.cell[0]**2 + xyz.cell[1]**2 + xyz.cell[2]**2)
        b = np.sqrt(xyz.cell[3]**2 + xyz.cell[4]**2 + xyz.cell[5]**2)
        c = np.sqrt(xyz.cell[6]**2 + xyz.cell[7]**2 + xyz.cell[8]**2)
        for atom in xyz.atoms:
            positions.append([atom.x / a, atom.y / b, atom.z / c])
            numbers.append(xyz.specie_labels[atom.name])
        structure = (lattice, positions, numbers)
        self.kpoints_seekpath = seekpath.get_path(structure)


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

    def to_input(self, fout):
        fout.write("\t\t\t&ELF_CUBE ON\n")
        fout.write("\t\t\t\tSTRIDE 1 1 1\n")
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
        self.params["FILENAME"] = "result.cube"
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
            self.active_space.to_print(fout)

        if self.band_structure.status == True:
            self.band_structure.to_print(fout)
        
        if self.efield_cube.status == True:
            self.efield_cube.to_print(fout)

        if self.electric_field_gradient.status == True:
            self.electric_field_gradient.to_print(fout)

        if self.elf_cube.status == True:
            self.elf_cube.to_print(fout)

        if self.e_density_cube.status == True:
            self.e_density_cube.to_print(fout)

        if self.external_potential_cube.status == True:
            self.external_potential_cube.to_print(fout)

        if self.pdos.status == True:
            self.pdos.to_print(fout)

        if self.lowdin.status == True:
            self.lowdin.to_print(fout)

        if self.mo.status == True:
            self.mo.to_print(fout)

        if self.moments.status == True:
            self.moments.to_print(fout)

        if self.mo_cubes.status == True:
            self.mo_cubes.to_print(fout)

        if self.mulliken.status == True:
            self.mulliken.to_print(fout)

        if self.stm.status == True:
            self.stm.to_print(fout)

        if self.subcell.status == True:
            self.subcell.to_print(fout)
       
        if self.tot_density_cube.status == True:
            self.tot_density_cube.to_print(fout)

        if self.v_hartree_cube.status == True:
            self.v_hartree_cube.to_print(fout)

        if self.v_xc_cube.status == True:
            self.v_xc_cube.to_print(fout)

        if self.wannier90.status == True:
            self.wannier90.to_print(fout)

        if self.wfn_mix.status == True:
            self.wfn_mix.to_print(fout)

        if self.xray_diffraction_spectrum.status == True:
            self.xray_diffraction_spectrum.to_print(fout)

        fout.write("\t\t&END PRINT\n")


    def set_params(self, params):
        for item in params:
            if len(item.spit("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            elif item.split("-")[2] == "BAND_STRUCTURE":
                self.band_structure.set_params({item: params[item]})
