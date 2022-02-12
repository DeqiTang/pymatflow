"""
in control of the properties calculation related parameters
"""
import numpy as np
import sys
import os
import shutil

from pymatflow.siesta.group import SiestaVariableGroup

"""
Usage:
"""


class SiestaProperties(SiestaVariableGroup):
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
    def __init__(self):
        """
        """
        super().__init__()
        self.xyz = None
        self.options = []

        self.bandlines = None

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

    def set_xyz(self, xyz):
        """
        :param xyz is the instance of base_xyz, passed from the calculation controller like static_run,
            it is used by seekpath to generate BandLines
        """
        self.xyz = xyz

    def to_string(self):
        out = ""
        out += super().to_string()
        
        out += "# ====================================\n"
        out += "# properties output related parameters\n"
        out += "# ====================================\n"
        out += "\n"
        #
        if 1 in self.options:
            out += self.to_string_pdos()
        if 2 in self.options:
            out += self.to_string_ldos()
        if 3 in self.options:
            out += self.to_string_band()
        if 4 in self.options:
            out += self.to_string_charge()
        if 5 in self.options:
            out += self.to_string_chem()
        if 6 in self.options:
            out += self.to_string_macro_polarization()
        if 7 in self.options:
            out += self.to_string_net_charge_dipole_elec_field()
        if 8 in self.options:
            out += self.to_string_optical()
        if 9 in self.options:
            out += self.to_string_wannier90()
        return out

    def to_string_pdos(self, emin=-20.0, emax=20.0, peak_width=0.5, npoint=2000, 
            kgrid_mp=[8, 0, 0, 0.5, 0, 8, 0, 0.5, 0, 0, 8, 0.5]):
        # 输出pdos时会同时输出total dos
        out = ""
        out += "%block ProjectedDensityOfStates\n"
        out += "%f %f %f %d eV\n" % (emin, emax, peak_width, npoint)
        out += "%endblock ProjectedDensityOfStates\n"
        out += "%block PDOS.kgrid.MonkhorstPack\n"
        out += "%d %d %d %.1f\n" % (kgrid_mp[0], kgrid_mp[1], kgrid_mp[2], kgrid_mp[3])
        out += "%d %d %d %.1f\n" % (kgrid_mp[4], kgrid_mp[5], kgrid_mp[6], kgrid_mp[7])
        out += "%d %d %d %.1f\n" % (kgrid_mp[8], kgrid_mp[9], kgrid_mp[10], kgrid_mp[11])
        out += "%endblock PDOS.kgrid.MonkhorstPack\n"
        out += "\n"
        return out

    def to_string_ldos(self, emin=-3.5, emax=0.0):
        out = ""
        out += "%block LocalDensityOfStates\n"
        out += "%f %f eV\n" % (emin, emax)
        out += "%endblock LocalDensityOfStates\n"
        return out

    def to_string_band(self, option=0):
        """
        option:
            0: using BandLines to specify kpoints for band calculation
            1: using BandPoints to specify kpoints for band calculation [not implemented now]

            self.bandlines:  please use ReciprocalLatticeVectors (Crystal) coordinates
                the high symmetry k point path used in bands structure calculation
                in format like this:
            
                [[kx, ky, kz, label, connect_indicator], ...] like [[0.0, 0.0, 0.0, 'GAMMA', 15], ...]
            
            if connect_indicator in a kpoint is an integer, then it will connect to the following point
            through the number of kpoints defined by connect_indicator.
            
            if connect_indicator in a kpoint is '|', then it will not connect to the following point,
        
        """
        out = ""
        #out += "BandLinesScale pi/a\n")
        out += "BandLinesScale ReciprocalLatticeVectors\n"
        out += "WriteKbands false\n" # do not write kbands to the standard out
        out += "WriteBands false\n" # do not write bands to the standard out
        if option == 0:
            out += "%block BandLines\n"
            out += "%d %f %f %f %s\n" % (1, self.bandlines[0][0], self.bandlines[0][1], self.bandlines[0][2], self.bandlines[0][3])
            for i in range(1, len(self.bandlines)):
                out += "%d %f %f %f %s\n" % (
                    self.bandlines[i-1][4] if type(self.bandlines[i-1][4]) == int else 1, # I want to set it to 0 like in qe, but it won't work normally, so set to 1 now.
                    self.bandlines[i][0], self.bandlines[i][1], self.bandlines[i][2], self.bandlines[i][3],
                    )
            out += "%endblock BandLines\n"
            out += "\n"
        elif option == 1:
            pass
        return out

    def to_string_charge(self):
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
        out = ""
        out += "Write.Denchar true\n"
        out += "SaveRho true\n"
        out += "SaveDeltaRho true\n"
        out += "SaveRhoXC true\n"
        out += "SaveElectricalstaticPotential true\n"
        out += "SaveNeutralAtomPotential true\n"
        out += "SaveTotalPotential true\n"
        out += "SaveIonicCharge true\n"
        out += "SaveTotalCharge true\n"
        # generated SystemLabel.BADER can be converted to cube format 
        # by Util/Grid/grid2cube(failed in my machine) or Util/Grid/g2c_ng
        # use of g2c_ng tool:
        # g2c_ng -g SystemLabel.BADER -s SystemLabel.STRUCT_OUT
        # and the get cube file can be further analyzed using
        # bader program: http://theory.cm.utexas.edu/bader/
        # bader SystemLabel.RHO.cube -ref SystemLabel.BADER.cube
        out += "SaveBaderCharge true\n"
        out += "AnalyzeChargeDensityOnly false\n"
        out += "SaveInitialChargeDensity false\n"
        out += "\n"
        return out

    def to_string_chem(self):
        # Mulliken charges and overlap populations
        out = ""
        out += "WriteMullikenPop 3\n"
        out += "MullikenInSCF false\n"
        out += "SpinInSCF true\n"
        # Voronoi and Hirshfeld atomic population analysis
        out += "WriteHirshfeldPop true\n"
        out += "WriteVoronoiPop     true\n"
        out += "PartialChargesAtEveryGeometry false\n"
        out += "PartialChargesAtEverySCFStep false\n"
        # Crystal-Orbital overlap and hamilton populations(COOP/COHP)
        out += "COOP.Write true\n"
        out += "#WFS.Energy.Min # default is minus infinity\n"   # default −∞
        out += "#WFS.Energy.Max # default is positive infinity\n"  # default ∞
        out += "\n"
        return out

    def to_string_macro_polarization(self):
        # Macroscopic Polarzation: using the geometric Berry phase approach
        # Refer to the official siesta manual about the meaning of the components
        # inside PolarizationGrids.
        out = ""
        out += "%block PolarizationGrids\n"
        for item in self.polarization_grids:
            out += "  %d %d %d %s\n" % (int(item.split()[0]), int(item.split()[1]), int(item.split()[2]), item.split()[3])
        out += "%endblock PolarizationGrids\n"
        out += "BornCharge true\n"
        out += "\n"
        return out

    def to_string_net_charge_dipole_elec_field(self):
        # Systems with net charge or dipole, and electric fields
        out = ""
        out += "NetCharge 0\n"
        out += "SimulateDoping true\n"
        out += "%block ExternalElectricField\n"
        out += "  %.3f %.3f %.3f V/Ang\n" % (float(self.external_electric_field[0]), float(self.external_electric_field[1]), float(self.external_electric_field[2]))
        out += "%endblock ExternalELectricField\n"
        out += "SlabDipoleCorrection true\n"
        # out += "%block Geometry.Hartree\n"
        # out += "%endblock Geometry.Hartree\n"
        # out += "%block Geometry.Charge\n"
        # out += "%endblock Geometry.Charge\n"
        out += "\n"
        return out

    def to_string_optical(self):
        # Optical properties
        out = ""
        out += "OpticalCalculation true\n"
        out += "Optical.Energy.Minimum %.3f Ry\n" % self.optical_energy_minimum
        out += "Optical.Energy.Maximum %.3f Ry\n" % self.optical_energy_maximum
        out += "Optical.Broaden %.3f Ry\n" % self.optical_broaden
        out += "Optical.Scissor %.3f Ry\n" % self.optical_scissor
        out += "#Optical.NumberOfBands # default: all bands\n" # default: all bands
        out += "%block Optical.Mesh\n"
        out += "  %d %d %d\n" % (self.optical_mesh[0], self.optical_mesh[1], self.optical_mesh[2])
        out += "%endblock Optical.Mesh\n"
        out += "Optical.OffsetMesh true\n"
        out += "Optical.PolarizationType %s\n" % self.optical_polarization_type # polycrystal, polarized, unpolarized
        out += "%block Optical.Vector\n"
        out += "  %.3f %.3f %.3f\n" % (self.optical_vector[0], self.optical_vector[1], self.optical_vector[2])
        out += "%endblock Optical.Vector\n"
        out += "\n"
        return out
    
    def to_string_wannier90(self):
        # wannier90: Maximally Localized Wannier Functions
        out = ""
        out += "Siesta2Wannier90.WriteMmm true\n"
        out += "Siesta2Wannier90.WriteAmm true\n"
        out += "Siesta2Wannier90.WriteEig true\n"
        out += "Siesta2Wannier90.WriteUnk true\n"
        out += "Siesta2Wannier90.UnkGrid1 %d\n" % self.wannier90_unkgrid[0]
        out += "Siesta2Wannier90.UnkGrid2 %d\n" % self.wannier90_unkgrid[1]
        out += "Siesta2Wannier90.UnkGrid3 %d\n" % self.wannier90_unkgrid[2]
        out += "Siesta2Wannier90.UnkGridBinary true\n"
        out += "#Siesta2Wannier90.NumberOfBands # default: occupied bands\n" # default: occupied bands
        out += "#Siesta2Wannier90.NumberOfBandsUp # default: Siesta2Wannier90.NumberOfBands\n"
        out += "#Siesta2Wannier90.NumberOfBandsDown # default: Siesta2Wannier90.NumberOfBands\n"
        out += "\n"

        return out

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

