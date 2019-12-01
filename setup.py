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
    name = "pymatflow",
    version = '0.1',
    keywords = ("Emulation, DFT"),
    description = "an emulation assistant",
    license = "",
    url = "",
    author_email = "emuxlab@163.com",
    packages = find_packages(),
    #data_files = [('noteseton/uis', ['noteseton/uis/app.ui', 'noteseton/uis/menubar.ui', 'noteseton/uis/appmenu.ui'])],
    include_package_data = True,
    platforms = "any",
    install_requires = ["pymatgen"],

    scripts = [
        "pymatflow/abinit/scripts/bands_abinit.py",
        "pymatflow/abinit/scripts/converge_ecut_abinit.py",
        "pymatflow/abinit/scripts/dos_abinit.py",
        "pymatflow/abinit/scripts/geo_opt_abinit.py",
        "pymatflow/abinit/scripts/md_abinit.py",
        "pymatflow/abinit/scripts/neb_abinit.py",
        "pymatflow/abinit/scripts/phonon_abinit.py",
        "pymatflow/abinit/scripts/phonon_with_phonopy_abinit.py",
        "pymatflow/abinit/scripts/abinit-md.py",
        "pymatflow/abinit/scripts/abinit-opt.py",
        "pymatflow/abinit/scripts/abinit-converge-ecut.py",
        "pymatflow/abinit/scripts/abinit-phonopy.py",
        "pymatflow/abinit/scripts/abinit-scf.py",
        "pymatflow/abinit/scripts/abinit-nscf.py",
        "pymatflow/abinit/scripts/abinit-neb.py",
        "pymatflow/abinit/scripts/abinit-dfpt.py",
        "pymatflow/abinit/scripts/abinit-band.py",
        "pymatflow/abinit/post/scripts/post-abinit-neb.py",
        "pymatflow/abinit/post/scripts/post-abinit-opt.py",
        "pymatflow/abinit/post/scripts/post-abinit-scf.py",
        "pymatflow/abinit/post/scripts/post-abinit-md.py",
        "pymatflow/cp2k/scripts/cp2k-aimd.py",
        "pymatflow/cp2k/scripts/bands_cp2k.py",
        "pymatflow/cp2k/scripts/converge_cutoff_cp2k.py",
        "pymatflow/cp2k/scripts/cp2k-converge-cutoff.py",
        "pymatflow/cp2k/scripts/converge_rel_cutoff_cp2k.py",
        "pymatflow/cp2k/scripts/cp2k-converge-rel-cutoff.py",
        "pymatflow/cp2k/scripts/electron_density_cp2k.py",
        "pymatflow/cp2k/scripts/geo_opt_cp2k.py",
        "pymatflow/cp2k/scripts/mp2_cp2k.py",
        "pymatflow/cp2k/scripts/pdos_cp2k.py",
        "pymatflow/cp2k/scripts/phonon_cp2k.py",
        "pymatflow/cp2k/scripts/single-point-cp2k.py",
        "pymatflow/cp2k/scripts/cp2k-geo-opt.py",
        "pymatflow/cp2k/scripts/cp2k-cell-opt.py",
        "pymatflow/cp2k/scripts/cp2k-scf.py",
        "pymatflow/cp2k/scripts/cp2k-scf-restart.py",
        "pymatflow/cp2k/scripts/cp2k-vib.py",
        "pymatflow/cp2k/scripts/cp2k-neb.py",
        "pymatflow/cp2k/scripts/cp2k-lr.py",
        "pymatflow/cp2k/post/scripts/post-single-point-cp2k.py",
        "pymatflow/cp2k/post/scripts/post-geo-opt-cp2k.py",
        "pymatflow/cp2k/post/scripts/post-pdos-without-convolute.py",
        "pymatflow/cp2k/post/scripts/post-cp2k-geo-opt.py",
        "pymatflow/cp2k/post/scripts/post-cp2k-cell-opt.py",
        "pymatflow/cp2k/post/scripts/post-cp2k-scf.py",
        "pymatflow/cp2k/post/scripts/post-cp2k-pdos.py",
        "pymatflow/cp2k/post/scripts/post-cp2k-bands.py",
        "pymatflow/cp2k/post/scripts/post-cp2k-converge.py",
        "pymatflow/cp2k/post/scripts/post-cp2k-md.py",
        "pymatflow/cp2k/post/scripts/cp2k_bs2csv.py",
        "pymatflow/dalton/scripts/dalton-single-point.py",
        "pymatflow/orca/single_point_orca.py",
        "pymatflow/orca/geo_opt_orca.py",
        "pymatflow/orca/converge_multiplicity_orca.py",
        "pymatflow/qe/scripts/qe-md.py",
        "pymatflow/qe/scripts/qe-vc-md.py",
        "pymatflow/qe/scripts/qe-converge-ecutwfc.py",
        "pymatflow/qe/scripts/qe-converge-ecutrho.py",
        "pymatflow/qe/scripts/qe-converge-kpoints.py",
        "pymatflow/qe/scripts/qe-converge-degauss.py",
        "pymatflow/qe/scripts/qe-scf-nscf.py",
        "pymatflow/qe/scripts/qe-scf-nscf-dos.py",
        "pymatflow/qe/scripts/qe-scf-nscf-bands.py",
        "pymatflow/qe/scripts/qe-nscf.py",
        "pymatflow/qe/scripts/qe-scf.py",
        "pymatflow/qe/scripts/qe-dos.py",
        "pymatflow/qe/scripts/qe-bands.py",
        "pymatflow/qe/scripts/qe-pdos.py",
        "pymatflow/qe/scripts/qe-epsilon.py",
        "pymatflow/qe/scripts/qe-ir-raman.py",
        "pymatflow/qe/scripts/qe-fermi-surface.py",
        "pymatflow/qe/scripts/qe-relax.py",
        "pymatflow/qe/scripts/qe-vc-relax.py",
        "pymatflow/qe/scripts/qe-neb.py",
        "pymatflow/qe/scripts/qe-turbo-davidson-spectrum.py",
        "pymatflow/qe/scripts/qe-turbo-lanczos-spectrum.py",
        "pymatflow/qe/scripts/qe-phx.py",
        "pymatflow/qe/scripts/qe-q2r.py",
        "pymatflow/qe/scripts/qe-matdyn.py",
        "pymatflow/qe/scripts/qe-dynmat.py",
        "pymatflow/qe/scripts/qe-plotband.py",
        "pymatflow/qe/scripts/qe-pp.py",
        "pymatflow/qe/scripts/qe-phonopy.py",
        "pymatflow/qe/post/scripts/post-qe-dos.py",
        "pymatflow/qe/post/scripts/post-qe-pdos.py",
        "pymatflow/qe/post/scripts/post-qe-neb.py",
        "pymatflow/qe/post/scripts/post-qe-relax.py",
        "pymatflow/qe/post/scripts/post-qe-vc-relax.py",
        "pymatflow/qe/post/scripts/post-qe-scf.py",
        "pymatflow/qe/post/scripts/post-qe-converge.py",
        "pymatflow/qe/post/scripts/post-qe-bands.py",
        "pymatflow/siesta/scripts/band_new_scf_siesta.py",
        "pymatflow/siesta/scripts/charge-density_potentials_siesta.py",
        "pymatflow/siesta/scripts/chem_analysis_siesta.py",
        "pymatflow/siesta/scripts/converge_ecut_siesta.py",
        "pymatflow/siesta/scripts/coord_from_xyz_to_siesta.py",
        "pymatflow/siesta/scripts/dos-pdos_new_scf_siesta.py",
        "pymatflow/siesta/scripts/geo_opt_siesta.py",
        "pymatflow/siesta/scripts/siesta-opt.py",
        "pymatflow/siesta/scripts/ldos_new_scf_siesta.py",
        "pymatflow/siesta/scripts/macro_polarization_siesta.py",
        "pymatflow/siesta/scripts/md_siesta.py",
        "pymatflow/siesta/scripts/siesta-md.py",
        "pymatflow/siesta/scripts/net-charge_dipole_elec-field_siesta.py",
        "pymatflow/siesta/scripts/optical_siesta.py",
        "pymatflow/siesta/scripts/phonon_siesta.py",
        "pymatflow/siesta/scripts/phonon_with_phonopy_siesta.py",
        "pymatflow/siesta/scripts/wannier90_siesta.py",
        "pymatflow/siesta/scripts/xyz2fdf.py",
        "pymatflow/siesta/scripts/siesta-converge-cutoff.py",
        "pymatflow/siesta/scripts/siesta-scf.py",
        "pymatflow/siesta/scripts/siesta-scf-restart.py",
        "pymatflow/siesta/scripts/siesta-ts.py",
        "pymatflow/siesta/post/scripts/post-siesta-pdos.py",
        "pymatflow/siesta/post/scripts/post-siesta-bands.py",
        "pymatflow/siesta/post/scripts/post-siesta-converge.py",
        "pymatflow/siesta/post/scripts/post-siesta-opt.py",
        "pymatflow/siesta/post/scripts/post-siesta-md.py",
        "pymatflow/tools/build_supercell.py",
        "pymatflow/tools/cluster_sphere.py",
        "pymatflow/tools/qe-fix-atoms.py",
        "pymatflow/tools/xyzinfo.py",
        "pymatflow/remote/scripts/thq.py",
        "pymatflow/remote/scripts/thcancel.py",
        "pymatflow/remote/scripts/thpull.py",
        "pymatflow/remote/scripts/thcmd.py",
        "pymatflow/remote/scripts/threport.py",
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



