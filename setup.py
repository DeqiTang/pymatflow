#!/usr/bin/env python
# _*_ coding: utf-8 _*_

from setuptools import setup, find_packages
from setuptools import Command
import os

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

setup(
    name = "emuhelper",
    version = '0.1',
    keywords = ("Emulation, DFT"),
    description = "an emulation assistant",
    license = "",
    url = "",
    author_email = "emulab@163.com",
    packages = find_packages(),
    #data_files = [('noteseton/uis', ['noteseton/uis/app.ui', 'noteseton/uis/menubar.ui', 'noteseton/uis/appmenu.ui'])],
    include_package_data = True,
    platforms = "any",
    install_requires = ["pymatgen"],

    scripts = [
        "emuhelper/abinit/bands_abinit.py",
        "emuhelper/abinit/converge_ecut_abinit.py",
        "emuhelper/abinit/dos_abinit.py",
        "emuhelper/abinit/geo_opt_abinit.py",
        "emuhelper/abinit/md_abinit.py",
        "emuhelper/abinit/neb_abinit.py",
        "emuhelper/abinit/phonon_abinit.py",
        "emuhelper/abinit/phonon_with_phonopy_abinit.py",
        "emuhelper/cp2k/scripts/aimd_cp2k.py",
        "emuhelper/cp2k/scripts/bands_cp2k.py",
        "emuhelper/cp2k/scripts/converge_cutoff_cp2k.py",
        "emuhelper/cp2k/scripts/converge_rel_cutoff_cp2k.py",
        "emuhelper/cp2k/scripts/electron_density_cp2k.py",
        "emuhelper/cp2k/scripts/geo_opt_cp2k.py",
        "emuhelper/cp2k/scripts/mp2_cp2k.py",
        "emuhelper/cp2k/scripts/pdos_cp2k.py",
        "emuhelper/cp2k/scripts/phonon_cp2k.py",
        "emuhelper/cp2k/scripts/single-point-cp2k.py",
        "emuhelper/cp2k/scripts/geo-opt-cp2k.py",
        "emuhelper/cp2k/post/scripts/post-single-point-cp2k.py",
        "emuhelper/cp2k/post/scripts/post-geo-opt-cp2k.py",
        "emuhelper/dalton/single_point_dalton.py",
        "emuhelper/orca/single_point_orca.py",
        "emuhelper/orca/geo_opt_orca.py",
        "emuhelper/orca/converge_multiplicity_orca.py",
        "emuhelper/qe/converge_ecutrho_pwscf.py",
        "emuhelper/qe/converge_ecutwfc_pwscf.py",
        "emuhelper/qe/geo_opt_not_vc_pwscf.py",
        "emuhelper/qe/geo_opt_vc_pwscf.py",
        "emuhelper/qe/phonon_with_phonopy_pwscf.py",
        "emuhelper/qe/phonon_with_phx_pwscf.py",
        "emuhelper/siesta/band_new_scf_siesta.py",
        "emuhelper/siesta/charge-density_potentials_siesta.py",
        "emuhelper/siesta/chem_analysis_siesta.py",
        "emuhelper/siesta/converge_ecut_siesta.py",
        "emuhelper/siesta/coord_from_xyz_to_siesta.py",
        "emuhelper/siesta/dos-pdos_new_scf_siesta.py",
        "emuhelper/siesta/geo_opt_siesta.py",
        "emuhelper/siesta/ldos_new_scf_siesta.py",
        "emuhelper/siesta/macro_polarization_siesta.py",
        "emuhelper/siesta/md_siesta.py",
        "emuhelper/siesta/net-charge_dipole_elec-field_siesta.py",
        "emuhelper/siesta/optical_siesta.py",
        "emuhelper/siesta/phonon_siesta.py",
        "emuhelper/siesta/phonon_with_phonopy_siesta.py",
        "emuhelper/siesta/wannier90_siesta.py",
        "emuhelper/siesta/xyz2fdf.py",
        "emuhelper/tools/build_supercell.py",
        "emuhelper/tools/cluster_sphere.py"
        ],
    #entry_points = {
    #    'console_scripts': [
    #        'noteseton = main:main'
    #    ]
    #},

    cmdclass = {
        'clean': CleanCommand,
    }
)



