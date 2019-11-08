#!/usr/bin/env python
# _*_ coding: utf-8 _*_


# ===============================
# CP2K / FORCE_EVAL / DFT / PRINT
# ===============================
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
        self.active_space = False
        self.adjmat_write = False
        self.band_structure = False
        self.basis_molopt_quantities = False
        self.basis_set_file = False
        self.derivatives = False
        self.dft_control_parameters = False
        self.efield_cube = False
        self.electric_field_gradient = False
        self.elf_cube = False
        self.energy_windows = False
        self.external_potential_cube = False
        self.e_density_cube = False
        self.gapw = False
        self.hirshfeld = False
        self.hyperfine_coupling_tensor = False
        self.implicit_psolver = False
        self.kinetic_energy = False
        self.kpoints = False
        self.ks_csr_write = False
        self.lowdin = False
        self.mao_analysis = False
        self.minbas_analysis = False
        self.mo = False
        self.moments = False
        self.mo_cubes = False
        self.mulliken = False
        self.neighbor_lists = False
        self.optimize_lri_basis = False
        self.overlap_condition = False
        self.pdos = False
        self.plus_u = False
        self.program_banner = False
        self.sccs = False
        self.stm = False
        self.subcell = False
        self.s_csr_write = False
        self.tot_density_cube = False
        self.v_hartree_cube = False
        self.v_xc_cube = False
        self.wannier90 = False
        self.wfn_mix = False
        self.xray_diffraction_spectrum = False

    def to_dft(self, fout):
        """
        fout: a file stream for writing
        """
        fout.write("\t\t&PRINT\n")
        for item in self.params:
            if self.params[item] is not None:
                fout.write("\t\t\t%s %s\n" % (item, self.params[item]))

        if self.band_structure == True:
            fout.write("\t\t\t&BAND_STRUCTURE\n")
            fout.write("\t\t\t\tADDED_MOS 10\n") 
            fout.write("\t\t\t\tFILE_NAME bands.bs\n")
            fout.write("\t\t\t\t&KPOINT_SET\n")
            fout.write("\t\t\t\t\tUNITS B_VECTOR\n")
            fout.write("\t\t\t\t\tSPECIAL_POINT GAMA 0.0000000000     0.0000000000     0.0000000000\n") # GAMA
            fout.write("\t\t\t\t\tSPECIAL_POINT X    0.5000000000     0.0000000000     0.5000000000\n") # X
            fout.write("\t\t\t\t\tSPECIAL_POINT W    0.5000000000     0.2500000000     0.7500000000\n") # W
            fout.write("\t\t\t\t\tSPECIAL_POINT GAMA 0.0000000000     0.0000000000     0.0000000000\n") # GAMA
            fout.write("\t\t\t\t\tSPECIAL_POINT U    0.6250000000     0.2500000000     0.6250000000\n") # U
            fout.write("\t\t\t\t\tSPECIAL_POINT W    0.5000000000     0.2500000000     0.7500000000\n") # W
            fout.write("\t\t\t\t\tNPOINTS 20\n")
            fout.write("\t\t\t\t&END KPOINT_SET\n")
            fout.write("\t\t\t&END BAND_STRUCTURE\n")
        
        if self.efield_cube == True:
            fout.write("\t\t\t&EFIELD_CUBE ON\n")
            fout.write("\t\t\t&END EFIELD_CUBE\n")

        if self.elf_cube == True:
            fout.write("\t\t\t&ELF_CUBE ON\n")
            fout.write("\t\t\t\tSTRIDE 1 1 1\n")
            fout.write("\t\t\t&END ELF_CUBE\n")

        if self.e_density_cube == True:
            fout.write("\t\t\t&E_DENSITY_CUBE ON\n")
            fout.write("\t\t\t&END E_DENSITY_CUBE\n")

        if self.pdos == True:
            fout.write("\t\t\t&PDOS\n")
            fout.write("\t\t\t\tNLUMO -1\n") # print all projected DOS available
            fout.write("\t\t\t\tCOMPONENTS\n") # split the density by quantum number
            fout.write("\t\t\t&END PDOS\n")

        if self.mo == True:
            fout.write("\t\t\t&MO\n")
            fout.write("\t\t\t&END MO\n")

        if self.moments == True:
            fout.write("\t\t\t&MOMENTS\n")
            fout.write("\t\t\t\tPERIODIC %s\n" % False)
            fout.write("\t\t\t&END MOMENTS\n")

        if self.mo_cubes == True:
            fout.write("\t\t\t&MO_CUBES\n")
            fout.write("\t\t\t\tNHOMO -1\n")
            fout.write("\t\t\t\tNLUMO -1\n")
            fout.write("\t\t\t&END MO_CUBES\n")
 
        if self.mulliken == True:
            fout.write("\t\t\t&MULLIKEN\n")
            fout.write("\t\t\t&END MULLIKEN\n")

        if self.stm == True:
            fout.write("\t\t\t&STM\n")
            fout.write("\t\t\t&END STM\n")
       
        if self.tot_density_cube == True:
            fout.write("\t\t\t&TOT_DENSITY_CUBE\n")
            fout.write("\t\t\t&END TOT_DENSITY_CUBE\n")

        if self.v_hartree_cube == True:
            fout.write("\t\t\t&V_HARTREE_CUBE\n")
            fout.write("\t\t\t&END V_HARTREE_CUBE\n")

        if self.v_xc_cube == True:
            fout.write("\t\t\t&V_XC_CUBE\n")
            fout.write("\t\t\t&END V_XC_CUBE\n")
        
        if self.xray_diffraction_spectrum == True:
            fout.write("\t\t\t&XRAY_DIFFRACTION_SPECTRUM\n")
            fout.write("\t\t\t&END XRAY_DIFFRACTION_SPECTRUM\n")

        fout.write("\t\t&END PRINT\n")

    def print_electron_density(self):
        self.electron_density = True

    def print_bands(self):
        self.bands = True
    def print_pdos(self):
        self.pdos = True
    def print_moments(self):
        self.moments = True

    def set_params(self, params):
        for item in params:
            if len(item.spit("-")) == 3:
                self.params[item.split("-")[-1]] = params[item]
            else:
                pass
