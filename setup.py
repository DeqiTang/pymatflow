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
    version = '0.0.2-beta',
    keywords = ("Ab intio ,DFT, workflow, input generation"),
    description = "An emulation assistant, input generation and manage for DFT programs",
    license = "MIT",
    url = "https://gitlab.com/deqitang/pymatflow",
    author_email = "pymatflow@163.com",
    packages = find_packages(),
    #data_files
    include_package_data = True,
    platforms = "any",
    install_requires = [],

    scripts = [
        "pymatflow/abinit/scripts/abinit-md.py",
        "pymatflow/abinit/scripts/abinit-opt.py",
        "pymatflow/abinit/scripts/abinit-converge-ecut.py",
        "pymatflow/abinit/scripts/abinit-phonopy.py",
        "pymatflow/abinit/scripts/abinit-scf.py",
        "pymatflow/abinit/scripts/abinit-nscf.py",
        "pymatflow/abinit/scripts/abinit-neb.py",
        "pymatflow/abinit/scripts/abinit-bands.py",
        "pymatflow/abinit/scripts/abinit-opt-cubic.py",
        "pymatflow/abinit/scripts/abinit-opt-tetragonal.py",
        "pymatflow/abinit/scripts/abinit-opt-hexagonal.py",
        "pymatflow/abinit/scripts/abinit-dfpt-elastic-piezo-dielec.py",
        "pymatflow/abinit/scripts/abinit-dfpt-phonon.py",
        "pymatflow/abinit/scripts/abinit-scf-nscf-dos-bands.py",
        "pymatflow/abinit/post/scripts/post-abinit-neb.py",
        "pymatflow/abinit/post/scripts/post-abinit-opt.py",
        "pymatflow/abinit/post/scripts/post-abinit-scf.py",
        "pymatflow/abinit/post/scripts/post-abinit-md.py",
        "pymatflow/abinit/post/scripts/post-abinit-phonopy.py",
        "pymatflow/abinit/post/scripts/post-abinit-dfpt-elastic-piezo-dielec.py",
        "pymatflow/cp2k/scripts/cp2k-aimd.py",
        "pymatflow/cp2k/scripts/cp2k-converge-cutoff.py",
        "pymatflow/cp2k/scripts/cp2k-converge-rel-cutoff.py",
        "pymatflow/cp2k/scripts/cp2k-converge-kpoints-manual.py",
        "pymatflow/cp2k/scripts/cp2k-converge-kpoints-auto.py",
        "pymatflow/cp2k/scripts/cp2k-geo-opt.py",
        "pymatflow/cp2k/scripts/cp2k-cell-opt.py",
        "pymatflow/cp2k/scripts/cp2k-scf.py",
        "pymatflow/cp2k/scripts/cp2k-mp2.py",
        "pymatflow/cp2k/scripts/cp2k-vib.py",
        "pymatflow/cp2k/scripts/cp2k-neb.py",
        "pymatflow/cp2k/scripts/cp2k-lr.py",
        "pymatflow/cp2k/scripts/cp2k-phonopy.py",
        "pymatflow/cp2k/scripts/cp2k-md-vib.py",
        "pymatflow/cp2k/scripts/cp2k-geo-opt-hexagonal.py",
        "pymatflow/cp2k/scripts/cp2k-geo-opt-tetragonal.py",
        "pymatflow/cp2k/scripts/cp2k-geo-opt-cubic.py",
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
        "pymatflow/cp2k/post/scripts/post-cp2k-phonopy.py",
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
        "pymatflow/qe/scripts/qe-plotband-for-matdyn.py",
        "pymatflow/qe/scripts/qe-pp.py",
        "pymatflow/qe/scripts/qe-phonopy.py",
        "pymatflow/qe/scripts/qe-molecularpdos.py",
        "pymatflow/qe/scripts/qe-pw-dielectric.py",
        "pymatflow/qe/scripts/qe-relax-tetragonal.py",
        "pymatflow/qe/scripts/qe-relax-cubic.py",
        "pymatflow/qe/scripts/qe-relax-hexagonal.py",
        "pymatflow/qe/scripts/qe-get-matdyn-qpoints-from-bands-calc.py",
        "pymatflow/qe/post/scripts/post-qe-dos.py",
        "pymatflow/qe/post/scripts/post-qe-pdos.py",
        "pymatflow/qe/post/scripts/post-qe-neb.py",
        "pymatflow/qe/post/scripts/post-qe-relax.py",
        "pymatflow/qe/post/scripts/post-qe-vc-relax.py",
        "pymatflow/qe/post/scripts/post-qe-scf.py",
        "pymatflow/qe/post/scripts/post-qe-converge.py",
        "pymatflow/qe/post/scripts/post-qe-bands.py",
        "pymatflow/qe/post/scripts/post-qe-matdyn.py",
        "pymatflow/qe/post/scripts/post-qe-phonopy.py",
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
        "pymatflow/siesta/scripts/siesta-phonopy.py",
        "pymatflow/siesta/scripts/wannier90_siesta.py",
        "pymatflow/siesta/scripts/xyz2fdf.py",
        "pymatflow/siesta/scripts/siesta-converge-cutoff.py",
        "pymatflow/siesta/scripts/siesta-scf.py",
        "pymatflow/siesta/scripts/siesta-scf-restart.py",
        "pymatflow/siesta/scripts/siesta-ts.py",
        "pymatflow/siesta/scripts/siesta-phonon.py",
        "pymatflow/siesta/scripts/siesta-opt-cubic.py",
        "pymatflow/siesta/scripts/siesta-opt-tetragonal.py",
        "pymatflow/siesta/scripts/siesta-opt-hexagonal.py",
        "pymatflow/siesta/post/scripts/post-siesta-pdos.py",
        "pymatflow/siesta/post/scripts/post-siesta-bands.py",
        "pymatflow/siesta/post/scripts/post-siesta-converge.py",
        "pymatflow/siesta/post/scripts/post-siesta-opt.py",
        "pymatflow/siesta/post/scripts/post-siesta-md.py",
        "pymatflow/siesta/post/scripts/post-siesta-phonopy.py",
        "pymatflow/scripts/xyz-build-supercell.py",
        "pymatflow/scripts/cluster_sphere.py",
        "pymatflow/scripts/xyz-fix-atoms.py",
        "pymatflow/scripts/xyzinfo.py",
        "pymatflow/scripts/cif-to-xyz-modified.py",
        "pymatflow/scripts/cif-to-pdb.py",
        "pymatflow/scripts/pdb-to-cif.py",
        "pymatflow/scripts/xyz-modified-to-crystal.py",
        "pymatflow/scripts/kpath-xyz-seekpath.py",
        "pymatflow/scripts/pot-from-xyz-modified.py",
        "pymatflow/remote/scripts/thq.py",
        "pymatflow/remote/scripts/thcancel.py",
        "pymatflow/remote/scripts/thpull.py",
        "pymatflow/remote/scripts/thcmd.py",
        "pymatflow/remote/scripts/sz-pull.py",
        "pymatflow/remote/scripts/sz-cmd.py",
        "pymatflow/remote/scripts/sz-del.py",
        "pymatflow/remote/scripts/sz-q.py",
        "pymatflow/remote/scripts/threport.py",
        "pymatflow/flow/scripts/flow-run.py",
        "pymatflow/flow/scripts/flow-pes-static-qe.py",
        "pymatflow/flow/scripts/flow-pes-relax-qe.py",
        "pymatflow/flow/scripts/flow-pes-relax-cp2k.py",
        ],
    entry_points = {
        'console_scripts': [
            'matflow = pymatflow.cmd.matflow:main',
            'mflow = pymatflow.cmd.matflow:main',
            'postflow = pymatflow.cmd.postflow:main',
            'pflow = pymatflow.cmd.postflow:main',
            'structflow = pymatflow.cmd.structflow:main',
            'sflow = pymatflow.cmd.structflow:main',
        ]
    },

    cmdclass = {
        'clean': CleanCommand,
    }
)
