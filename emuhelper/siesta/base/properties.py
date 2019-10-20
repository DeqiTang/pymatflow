#!/usr/bin/evn python
# _*_ coding: utf-8 _*_

import numpy as np
import sys
import os
import shutil
import pymatgen as mg



"""
Usage:
"""


class siesta_properties:
    """
    self.option:
        0: None
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
        self.params = {
                }
        self.option = 0 

    def to_fdf(self, fout):
        for item in self.params:
            if self.params[item] is not None:
                fout.write("%s %s\n" % (item, str(self.params[item])))
        #
        if self.option == 0:
            return
        elif self.option == 1:
            self.to_fdf_pdos(fout)
        elif self.option == 2:
            self.to_fdf_ldos(fout)
        elif self.option == 3:
            self.to_fdf_band(fout)
        elif self.option == 4:
            self.to_fdf_charge(fout)
        elif self.option == 5:
            self.to_fdf_chem(fout)
        elif self.option == 6:
            self.to_fdf_macro_polarization(fout)
        elif self.option == 7:
            self.to_fdf_net_charge_dipole_elec_field(fout)
        elif self.option == 8:
            self.to_fdf_optical(fout)
        elif self.option == 9:
            self.to_fdf_wannier90(fout)

    def to_fdf_pdos(self, fout, emin=-20.0, emax=20.0, peak_width=0.1, npoint=2000):
        # 输出pdos时会同时输出total dos
        fout.write("%block ProjectedDensityOfStates\n")
        fout.write("%f %f %f %d eV\n" % (emin, emax, peak_width, npoint))
        fout.write("%endblock ProjectedDensityOfStates\n")
        fout.write("%block PDOS.kgrid.MonkhorstPack\n")
        fout.write("10 0 0 0.5\n")
        fout.write("0 10 0 0.5\n")
        fout.write("0 0 10 0.5\n")
        fout.write("%endblock PDOS.kgrid.MonkhorstPack\n")
        fout.write("\n")

    def to_fdf_ldos(self, fout, emin=-3.5, emax=0.0):
        fout.write("%block LocalDensityOfStates\n")
        fout.write("%f %f eV\n" % (emin, emax))
        fout.write("%endblock LocalDensityOfStates\n")

    def to_fdf_band(self, fout):        
        fout.write("BandLinesScale pi/a\n")
        fout.write("%block BandLines\n")
        fout.write("1 1.0 1.0 1.0 L # begin at L\n")
        fout.write("20 0.0 0.0 0.0 \Gamma # 20 points from L to gamma\n")
        fout.write("25 2.0 0.0 0.0 X # 25 points from gamma to X \n")
        fout.write("30 2.0 2.0 2.0 \Gamma # 30 points from X to gamma\n")
        fout.write("%endblock BandLines\n")
        fout.write("WriteKbands false\n") # do not write kbands to the standard out
        fout.write("WriteBands false\n") # do not write bands to the standard out
        fout.write("\n")

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
        fout.write("WriteVoronoiPop true\n")
        fout.write("PartialChargesAtEveryGeometry false\n")
        fout.write("PartialChargesAtEverySCFStep false\n")
        # Crystal-Orbital overlap and hamilton populations(COOP/COHP)
        fout.write("COOP.Write true\n")
        fout.write("#WFS.Energy.Min")   # default −∞
        fout.write("#WFS.Energy.Max")   # default ∞
        fout.write("\n")

    def to_fdf_macro_polarization(self, fout):
        # Macroscopic Polarzation
        fout.write("%block PolarizationGrids\n")
        fout.write("  10 3 4 yes\n")
        fout.write("  2 20 2 no\n")
        fout.write("  4 4 15\n")
        fout.write("%endblock PolarizationGrids\n")
        fout.write("BornCharge true\n")
        fout.write("\n")

    def to_fdf_net_charge_dipole_elec_field(self, fout):
        # Systems with net charge or dipole, and electric fields
        fout.write("NetCharge 1.0\n")
        fout.write("SimulateDoping true\n")
        fout.write("%block ExternalElectricField\n")
        fout.write("  0.000 0.000 0.50 V/Ang\n")
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
        fout.write("Optical.Energy.Minimum 0 Ry\n")
        fout.write("Optical.Energy.Maximum 10 Ry\n")
        fout.write("Optical.Broaden 0 Ry\n")
        fout.write("Optical.Scissor 0 Ry\n")
        fout.write("#Optical.NumberOfBands # default: all bands\n") # default: all bands
        fout.write("%block Optical.Mesh\n")
        fout.write("  5 5 5\n")
        fout.write("%endblock Optical.Mesh\n")
        fout.write("Optical.OffsetMesh true\n")
        fout.write("Optical.PolarizationType unpolarized\n") # polycrystal, polarized, unpolarized
        fout.write("%block Optical.Vector\n")
        fout.write("  1.0 1.0 0.0\n")
        fout.write("%endblock Optical.Vector\n")
        fout.write("\n")
    
    def to_fdf_wannier90(self, fout):
        # wannier90: Maximally Localized Wannier Functions
        fout.write("Siesta2Wannier90.WriteMmm true\n")
        fout.write("Siesta2Wannier90.WriteAmm true\n")
        fout.write("Siesta2Wannier90.WriteEig true\n")
        fout.write("Siesta2Wannier90.WriteUnk true\n")
        fout.write("Siesta2Wannier90.UnkGrid1 10\n")
        fout.write("Siesta2Wannier90.UnkGrid2 10\n")
        fout.write("Siesta2Wannier90.UnkGrid3 10\n")
        fout.write("Siesta2Wannier90.UnkGridBinary true\n")
        fout.write("#Siesta2Wannier90.NumberOfBands # default: occupied bands\n") # default: occupied bands
        fout.write("#Siesta2Wannier90.NumberOfBandsUp # default: Siesta2Wannier90.NumberOfBands")
        fout.write("#Siesta2Wannier90.NumberOfBandsDown # default: Siesta2Wannier90.NumberOfBands")
        fout.write("\n")

    #

