#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg

import seekpath


"""
Usage:
"""


class siesta_properties:
    """
    self.option:
        1: PDOS
        2: LDOS
        3: Bands
        4: Charge Density
        5: Chemical analysis
        6: Macro Polarization
        7: Net Charge Dipole Electric Field
        8: Optical
        9: Wannier90


    """
    def __init__(self, xyz):
        """
        xyz is the instance of base_xyz, passed from the calculation controller like static_run,
        it is used by seekpath to generate BandLines
        """
        self.xyz = xyz
        self.params = {
                }
        self.options = []

        self.polarization_grids = [
                "10 3 3 no",
                "2 20 2 no",
                "4 4 15 no",
                ]
        self.external_electric_field = [0.0, 0.0, 0.5]
        self.optical_energy_minimum = 0 # Ry
        self.optical_energy_maximum = 10 # Ry
        self.optical_broaden = 0 # Ry
        self.optical_scissor = 0 # Ry
        self.optical_mesh = [5, 5, 5]
        self.optical_polarization_type = "unpolarized" # polarized, unpolarized, polycrystal
        self.optical_vector = [1.0, 0.0, 0.5]
        self.wannier90_unkgrid = [10, 10, 10]

    def to_fdf(self, fout):
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item, str(self.params[item])))
        
        fout.write("# ====================================\n")
        fout.write("# properties output related parameters\n")
        fout.write("# ====================================\n")
        fout.write("\n")
        #
        if 1 in self.options:
            self.to_fdf_pdos(fout)
        if 2 in self.options:
            self.to_fdf_ldos(fout)
        if 3 in self.options:
            self.to_fdf_band(fout)
        if 4 in self.options:
            self.to_fdf_charge(fout)
        if 5 in self.options:
            self.to_fdf_chem(fout)
        if 6 in self.options:
            self.to_fdf_macro_polarization(fout)
        if 7 in self.options:
            self.to_fdf_net_charge_dipole_elec_field(fout)
        if 8 in self.options:
            self.to_fdf_optical(fout)
        if 9 in self.options:
            self.to_fdf_wannier90(fout)

    def to_fdf_pdos(self, fout, emin=-20.0, emax=20.0, peak_width=0.5, npoint=2000, 
            kgrid_mp=[8, 0, 0, 0.5, 0, 8, 0, 0.5, 0, 0, 8, 0.5]):
        # 输出pdos时会同时输出total dos
        fout.write("%block ProjectedDensityOfStates\n")
        fout.write("%f %f %f %d eV\n" % (emin, emax, peak_width, npoint))
        fout.write("%endblock ProjectedDensityOfStates\n")
        fout.write("%block PDOS.kgrid.MonkhorstPack\n")
        fout.write("%d %d %d %.1f\n" % (kgrid_mp[0], kgrid_mp[1], kgrid_mp[2], kgrid_mp[3]))
        fout.write("%d %d %d %.1f\n" % (kgrid_mp[4], kgrid_mp[5], kgrid_mp[6], kgrid_mp[7]))
        fout.write("%d %d %d %.1f\n" % (kgrid_mp[8], kgrid_mp[9], kgrid_mp[10], kgrid_mp[11]))
        fout.write("%endblock PDOS.kgrid.MonkhorstPack\n")
        fout.write("\n")

    def to_fdf_ldos(self, fout, emin=-3.5, emax=0.0):
        fout.write("%block LocalDensityOfStates\n")
        fout.write("%f %f eV\n" % (emin, emax))
        fout.write("%endblock LocalDensityOfStates\n")

    def to_fdf_band(self, fout, option=0):
        """
        option:
            0: using BandLines to specify kpoints for band calculation
            1: using BandPoints to specify kpoints for band calculation [not implemented now]
        self.bandlines:
            a list of strings in format like below:
            ["1 0.000 0.000 0.000 \Gamma", ..., "20 1.0 1.0 1.0 L"]
        self. bandpoints:
            a list of strings in format like below:
            ["0.000 0.000 0.000", "1.000 0.000 0.000"]
        """
        fout.write("BandLinesScale pi/a\n")
        fout.write("WriteKbands false\n") # do not write kbands to the standard out
        fout.write("WriteBands false\n") # do not write bands to the standard out
        if option == 0:
            # --------------
            # using seekpath
            # --------------
            lattice = [self.xyz.cell[0:3], self.xyz.cell[3:6], self.xyz.cell[6:9]]
            positions = []
            numbers = []
            a = np.sqrt(self.xyz.cell[0]**2 + self.xyz.cell[1]**2 + self.xyz.cell[2]**2)
            b = np.sqrt(self.xyz.cell[3]**2 + self.xyz.cell[4]**2 + self.xyz.cell[5]**2)
            c = np.sqrt(self.xyz.cell[6]**2 + self.xyz.cell[7]**2 + self.xyz.cell[8]**2)
            for atom in self.xyz.atoms:
                positions.append([atom.x / a, atom.y / b, atom.z / c])
                numbers.append(self.xyz.specie_labels[atom.name])
            structure = (lattice, positions, numbers)
            self.kpoints_seekpath = seekpath.get_path(structure)

            fout.write("%block BandLines\n")
            point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][0][0]]
            fout.write("%d %f %f %f %s\n" % (1, point[0], point[1], point[2], self.kpoints_seekpath["path"][0][0]))
            point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][0][1]]
            fout.write("%d %f %f %f %s\n" % (20, point[0], point[1], point[2], self.kpoints_seekpath["path"][0][1]))
            for i in range(1, len(self.kpoints_seekpath["path"])):
                if self.kpoints_seekpath["path"][i][0] == self.kpoints_seekpath["path"][i-1][1]:
                    point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][1]]
                    fout.write("%d %f %f %f %s\n" % (20, point[0], point[1], point[2], self.kpoints_seekpath["path"][i][1]))
                else:
                    point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][0]]
                    fout.write("%d %f %f %f %s\n" % (20, point[0], point[1], point[2], self.kpoints_seekpath["path"][i][0]))
                    point = self.kpoints_seekpath["point_coords"][self.kpoints_seekpath["path"][i][1]]
                    fout.write("%d %f %f %f %s\n" % (20, point[0], point[1], point[2], self.kpoints_seekpath["path"][i][1]))
            fout.write("%endblock BandLines\n")
            fout.write("\n")
        elif option == 1:
            #fout.write("%block BandPoints\n")
            #for point in self.bandpoints:
            #    fout.write("%.3f %.3f %.3f\n" % (float(point.split()[0]), float(point.split()[1]), float(point.split()[2])))
            #fout.write("%endblock BandPoints\n")
            pass

    def to_fdf_charge(self, fout):
        # Output of charge densities and potentials on the grid
        """
        they can be processed by:
            Util/Grid
            Util/Contrib/APostnikov
            Util/Denchar
        Note:
            The program denchar in Util/Denchar can generate charge-density and wavefunction information
            in real space.
        """
        fout.write("Write.Denchar true\n")
        fout.write("SaveRho true\n")
        fout.write("SaveDeltaRho true\n")
        fout.write("SaveRhoXC true\n")
        fout.write("SaveElectricalstaticPotential true\n")
        fout.write("SaveNeutralAtomPotential true\n")
        fout.write("SaveTotalPotential true\n")
        fout.write("SaveIonicCharge true\n")
        fout.write("SaveTotalCharge true\n")
        # generated SystemLabel.BADER can be converted to cube format 
        # by Util/Grid/grid2cube(failed in my machine) or Util/Grid/g2c_ng
        # use of g2c_ng tool:
        # g2c_ng -g SystemLabel.BADER -s SystemLabel.STRUCT_OUT
        # and the get cube file can be further analyzed using
        # bader program: http://theory.cm.utexas.edu/bader/
        # bader SystemLabel.RHO.cube -ref SystemLabel.BADER.cube
        fout.write("SaveBaderCharge true\n") 
        fout.write("AnalyzeChargeDensityOnly false\n")
        fout.write("SaveInitialChargeDensity false\n")
        fout.write("\n")
    
    def to_fdf_chem(self, fout):
        # Mulliken charges and overlap populations
        fout.write("WriteMullikenPop 3\n")
        fout.write("MullikenInSCF false\n")
        fout.write("SpinInSCF true\n")
        # Voronoi and Hirshfeld atomic population analysis
        fout.write("WriteHirshfeldPop true\n")
        fout.write("WriteVoronoiPop     true\n")
        fout.write("PartialChargesAtEveryGeometry false\n")
        fout.write("PartialChargesAtEverySCFStep false\n")
        # Crystal-Orbital overlap and hamilton populations(COOP/COHP)
        fout.write("COOP.Write true\n")
        fout.write("#WFS.Energy.Min # default is minus infinity\n")   # default −∞
        fout.write("#WFS.Energy.Max # default is positive infinity\n")   # default ∞
        fout.write("\n")

    def to_fdf_macro_polarization(self, fout):
        # Macroscopic Polarzation: using the geometric Berry phase approach
        # Refer to the official siesta manual about the meaning of the components
        # inside PolarizationGrids.
        fout.write("%block PolarizationGrids\n")
        for item in self.polarization_grids:
            fout.write("  %d %d %d %s\n" % (int(item.split()[0]), int(item.split()[1]), int(item.split()[2]), item.split()[3]))
        fout.write("%endblock PolarizationGrids\n")
        fout.write("BornCharge true\n")
        fout.write("\n")

    def to_fdf_net_charge_dipole_elec_field(self, fout):
        # Systems with net charge or dipole, and electric fields
        fout.write("NetCharge 0\n")
        fout.write("SimulateDoping true\n")
        fout.write("%block ExternalElectricField\n")
        fout.write("  %.3f %.3f %.3f V/Ang\n" % (float(self.external_electric_field[0]), float(self.external_electric_field[1]), float(self.external_electric_field[2])))
        fout.write("%endblock ExternalELectricField\n")
        fout.write("SlabDipoleCorrection true\n")
        #fout.write("%block Geometry.Hartree\n")
        #fout.write("%endblock Geometry.Hartree\n")
        #fout.write("%block Geometry.Charge\n")
        #fout.write("%endblock Geometry.Charge\n")
        fout.write("\n")

    def to_fdf_optical(self, fout):
        # Optical properties
        fout.write("OpticalCalculation true\n")
        fout.write("Optical.Energy.Minimum %.3f Ry\n" % self.optical_energy_minimum)
        fout.write("Optical.Energy.Maximum %.3f Ry\n" % self.optical_energy_maximum)
        fout.write("Optical.Broaden %.3f Ry\n" % self.optical_broaden)
        fout.write("Optical.Scissor %.3f Ry\n" % self.optical_scissor)
        fout.write("#Optical.NumberOfBands # default: all bands\n") # default: all bands
        fout.write("%block Optical.Mesh\n")
        fout.write("  %d %d %d\n" % (self.optical_mesh[0], self.optical_mesh[1], self.optical_mesh[2]))
        fout.write("%endblock Optical.Mesh\n")
        fout.write("Optical.OffsetMesh true\n")
        fout.write("Optical.PolarizationType %s\n" % self.optical_polarization_type) # polycrystal, polarized, unpolarized
        fout.write("%block Optical.Vector\n")
        fout.write("  %.3f %.3f %.3f\n" % (self.optical_vector[0], self.optical_vector[1], self.optical_vector[2]))
        fout.write("%endblock Optical.Vector\n")
        fout.write("\n")
    
    def to_fdf_wannier90(self, fout):
        # wannier90: Maximally Localized Wannier Functions
        fout.write("Siesta2Wannier90.WriteMmm true\n")
        fout.write("Siesta2Wannier90.WriteAmm true\n")
        fout.write("Siesta2Wannier90.WriteEig true\n")
        fout.write("Siesta2Wannier90.WriteUnk true\n")
        fout.write("Siesta2Wannier90.UnkGrid1 %d\n" % self.wannier90_unkgrid[0])
        fout.write("Siesta2Wannier90.UnkGrid2 %d\n" % self.wannier90_unkgrid[1])
        fout.write("Siesta2Wannier90.UnkGrid3 %d\n" % self.wannier90_unkgrid[2])
        fout.write("Siesta2Wannier90.UnkGridBinary true\n")
        fout.write("#Siesta2Wannier90.NumberOfBands # default: occupied bands\n") # default: occupied bands
        fout.write("#Siesta2Wannier90.NumberOfBandsUp # default: Siesta2Wannier90.NumberOfBands\n")
        fout.write("#Siesta2Wannier90.NumberOfBandsDown # default: Siesta2Wannier90.NumberOfBands\n")
        fout.write("\n")

    #

    def set_params(self,
        #bandlines = [
        #        "1 0.000 0.000 0.000 \Gamma",
        #        "20 1.000 1.000 1.000 L",
        #        "20 2.000 0.000 0.000 X",
        #        ""
        #        ],
        #bandpoints = [
        #        "0.000 0.000 0.000",
        #        "1.000 0.000 0.000",
        #        "0.500 0.500 0.500",
        #        ],
        polarization_grids = [
                "10 3 3 no",
                "2 20 2 no",
                "4 4 15 no",
                ],
        external_electric_field = [0.0, 0.0, 0.5],
        optical_energy_minimum = 0, # Ry
        optical_energy_maximum = 10, # Ry
        optical_broaden = 0, # Ry
        optical_scissor = 0, # Ry
        optical_mesh = [5, 5, 5],
        optical_polarization_type = "unpolarized", # polarized, unpolarized, polycrystal
        optical_vector = [1.0, 0.0, 0.5],
        wannier90_unkgrid = [10, 10, 10]):
        """
        set params needed by various properties calculation
        """
        #self.bandlines = bandlines
        #self.bandpoints = bandpoints
        self.polarization_grids = polarization_grids
        self.external_electric_field = external_electric_field
        self.optical_energy_minimum = optical_energy_minimum
        self.optical_energy_maximum = optical_energy_maximum
        self.optical_broaden = optical_broaden
        self.optical_scissor = optical_scissor
        self.optical_mesh = optical_mesh
        self.optical_polarization_type = optical_polarization_type
        self.optical_vector = optical_vector
        self.wannier90_unkgrid = wannier90_unkgrid

