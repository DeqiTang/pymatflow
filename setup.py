#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import sys
import re
import subprocess
from glob import glob

import setuptools
#from setuptools import setup
from skbuild import setup  # replace setuptools setup
from setuptools import find_packages, Extension
from setuptools import Command
from setuptools.command.build_ext import build_ext

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')
        os.system("rm -rf _skbuild")


# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}

# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        # required for auto-detection & inclusion of auxiliary "native" libs
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={}".format(extdir),
            "-DPYTHON_EXECUTABLE={}".format(sys.executable),
            "-DCMAKE_BUILD_TYPE={}".format(cfg),  # not used on MSVC, but no harm
        ]
        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # In this example, we pass in the version to C++. You might not need to.
        cmake_args += [
            "-DEXAMPLE_VERSION_INFO={}".format(self.distribution.get_version())
        ]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator:
                try:
                    import ninja  # noqa: F401

                    cmake_args += ["-GNinja"]
                except ImportError:
                    pass

        else:

            # Single config generators are handled "normally"
            single_config = any(x in cmake_generator for x in {"NMake", "Ninja"})

            # CMake allows an arch-in-generator style for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
                ]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += ["-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call, not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += ["-j{}".format(self.parallel)]

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args, cwd=self.build_temp
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_temp
        )

# -----------------------------
# make sure pybind11 is working
# -----------------------------
try:
    import pybind11
except:
    subprocess.run(["pip3", "install", "--user", "pybind11[global]==2.7.1"])
os.environ["PATH"] = os.environ["PATH"] + ":%s" % os.path.join(os.path.expanduser("~"), ".local/bin")
sub = subprocess.run(["pybind11-config", "--includes"], stdout=subprocess.PIPE)
pybind11_include_paths = sub.stdout.decode().replace("\n", "")
# there might be more than one dir in output of pybind11-config --includes
list_paths = pybind11_include_paths.replace(" ", "").split("-I")
while "" in list_paths:
    list_paths.remove("")
if "CPLUS_INCLUDE_PATH" not in os.environ:
    os.environ["CPLUS_INCLUDE_PATH"] = ":".join(list_paths)
else:
    os.environ["CPLUS_INCLUDE_PATH"] = os.environ["CPLUS_INCLUDE_PATH"] + ":%s" % ":".join(list_paths)

# ---------------
# fortran support
# ---------------
#ar crv libbasic.a crystal_mod.o constant_mod.o cube_mod.o
# gfortran -shared -o libbasic.so crystal_mod.o constant_mod.o cube_mod.o
from numpy.distutils.core import Extension as npExtension
f_ext_modules = [
    npExtension(
        # module name:
        'cube2vtk',
        # source file: 
        sources = [
            'fortran/pyx/cube_to_vtk.pyx',
            'fortran/src/cube_to_vtk_c_binding_mod.f90'
        ],
        # other compile args for gcc
        #extra_compile_args=['-fPIC', '-O3'],
        # other files to link to
        #extra_link_args=['fortran/src/libbasic.so'],
        extra_link_args=['fortran/src/libbasic.so'],
        language="fortran"
    ),
]



setup(
    name = "pymatflow",
    version = '0.1.1a1',
    ## python3 setup.py build sdist bdist_wheel
    ## twine upload dist/*
    keywords = ("Ab intio ,DFT, workflow, input generation"),
    description = "An emulation assistant, input generation and manage for DFT programs",
    license = "MIT",
    url = "https://gitlab.com/deqitang/pymatflow",
    author_email = "pymatflow@163.com",
    packages = find_packages(),
    #data_files
    include_package_data = True,
    platforms = "any",
    python_requires = ">=3.0",
    install_requires = ["ase", "pybind11", "numpy", "scipy", "matplotlib"],
    # -----------------------------------------------
    # cpp extension using native setuptools (working)
    # -----------------------------------------------
    #ext_modules=[CMakeExtension(name="pyaskit", sourcedir="cpp")],
    #cmdclass = {
    #    'clean': CleanCommand,
    #    'build_ext': CMakeBuild,
    #},    
    # ------------------------------------------
    # cpp extension using scikit-build (working)
    # ------------------------------------------
    cmake_source_dir="cpp", # where CMakeLists.txt exists
    cmake_install_dir="pymatflow/cpp", # from pymatflow.cpp import pyaskit
    #cmake_install_dir="cpp/cppmatflow", # from cpp/cppmatflow import pyaskit
    cmdclass = {
        'clean': CleanCommand,
    },
    #cmake_args=["-DXXX=XX"],
    # -------
    # fortran
    # -------
    #ext_modules = f_ext_modules,
    # ----
    scripts = [
        #"pymatflow/abinit/scripts/abinit-md.py",
        #"pymatflow/abinit/scripts/abinit-opt.py",
        #"pymatflow/abinit/scripts/abinit-converge-ecut.py",
        #"pymatflow/abinit/scripts/abinit-phonopy.py",
        #"pymatflow/abinit/scripts/abinit-scf.py",
        #"pymatflow/abinit/scripts/abinit-nscf.py",
        #"pymatflow/abinit/scripts/abinit-neb.py",
        #"pymatflow/abinit/scripts/abinit-bands.py",
        #"pymatflow/abinit/scripts/abinit-opt-cubic.py",
        #"pymatflow/abinit/scripts/abinit-opt-tetragonal.py",
        #"pymatflow/abinit/scripts/abinit-opt-hexagonal.py",
        #"pymatflow/abinit/scripts/abinit-dfpt-elastic-piezo-dielec.py",
        #"pymatflow/abinit/scripts/abinit-dfpt-phonon.py",
        #"pymatflow/abinit/scripts/abinit-scf-nscf-dos-bands.py",
        #"pymatflow/abinit/post/scripts/post-abinit-neb.py",
        #"pymatflow/abinit/post/scripts/post-abinit-opt.py",
        #"pymatflow/abinit/post/scripts/post-abinit-scf.py",
        #"pymatflow/abinit/post/scripts/post-abinit-md.py",
        #"pymatflow/abinit/post/scripts/post-abinit-phonopy.py",
        #"pymatflow/abinit/post/scripts/post-abinit-dfpt-elastic-piezo-dielec.py",
        #"pymatflow/cp2k/scripts/cp2k-aimd.py",
        #"pymatflow/cp2k/scripts/cp2k-converge-cutoff.py",
        #"pymatflow/cp2k/scripts/cp2k-converge-rel-cutoff.py",
        #"pymatflow/cp2k/scripts/cp2k-converge-kpoints-manual.py",
        #"pymatflow/cp2k/scripts/cp2k-converge-kpoints-auto.py",
        #"pymatflow/cp2k/scripts/cp2k-geo-opt.py",
        #"pymatflow/cp2k/scripts/cp2k-cell-opt.py",
        #"pymatflow/cp2k/scripts/cp2k-scf.py",
        #"pymatflow/cp2k/scripts/cp2k-mp2.py",
        #"pymatflow/cp2k/scripts/cp2k-vib.py",
        #"pymatflow/cp2k/scripts/cp2k-neb.py",
        #"pymatflow/cp2k/scripts/cp2k-lr.py",
        #"pymatflow/cp2k/scripts/cp2k-phonopy.py",
        #"pymatflow/cp2k/scripts/cp2k-md-vib.py",
        #"pymatflow/cp2k/scripts/cp2k-geo-opt-hexagonal.py",
        #"pymatflow/cp2k/scripts/cp2k-geo-opt-tetragonal.py",
        #"pymatflow/cp2k/scripts/cp2k-geo-opt-cubic.py",
        #"pymatflow/cp2k/post/scripts/post-single-point-cp2k.py",
        #"pymatflow/cp2k/post/scripts/post-geo-opt-cp2k.py",
        #"pymatflow/cp2k/post/scripts/post-pdos-without-convolute.py",
        #"pymatflow/cp2k/post/scripts/post-cp2k-geo-opt.py",
        #"pymatflow/cp2k/post/scripts/post-cp2k-cell-opt.py",
        #"pymatflow/cp2k/post/scripts/post-cp2k-scf.py",
        #"pymatflow/cp2k/post/scripts/post-cp2k-pdos.py",
        #"pymatflow/cp2k/post/scripts/post-cp2k-bands.py",
        #"pymatflow/cp2k/post/scripts/post-cp2k-converge.py",
        #"pymatflow/cp2k/post/scripts/post-cp2k-md.py",
        #"pymatflow/cp2k/post/scripts/post-cp2k-phonopy.py",
        #"pymatflow/cp2k/post/scripts/cp2k_bs2csv.py",
        #"pymatflow/dalton/scripts/dalton-single-point.py",
        #"pymatflow/orca/single_point_orca.py",
        #"pymatflow/orca/geo_opt_orca.py",
        #"pymatflow/orca/converge_multiplicity_orca.py",
        #"pymatflow/qe/scripts/qe-md.py",
        #"pymatflow/qe/scripts/qe-vc-md.py",
        #"pymatflow/qe/scripts/qe-converge-ecutwfc.py",
        #"pymatflow/qe/scripts/qe-converge-ecutrho.py",
        #"pymatflow/qe/scripts/qe-converge-kpoints.py",
        #"pymatflow/qe/scripts/qe-converge-degauss.py",
        #"pymatflow/qe/scripts/qe-nscf.py",
        #"pymatflow/qe/scripts/qe-scf.py",
        #"pymatflow/qe/scripts/qe-dos.py",
        #"pymatflow/qe/scripts/qe-bands.py",
        #"pymatflow/qe/scripts/qe-pdos.py",
        #"pymatflow/qe/scripts/qe-epsilon.py",
        #"pymatflow/qe/scripts/qe-ir-raman.py",
        #"pymatflow/qe/scripts/qe-fermi-surface.py",
        #"pymatflow/qe/scripts/qe-relax.py",
        #"pymatflow/qe/scripts/qe-vc-relax.py",
        #"pymatflow/qe/scripts/qe-neb.py",
        #"pymatflow/qe/scripts/qe-turbo-davidson-spectrum.py",
        #"pymatflow/qe/scripts/qe-turbo-lanczos-spectrum.py",
        #"pymatflow/qe/scripts/qe-phx.py",
        #"pymatflow/qe/scripts/qe-q2r.py",
        #"pymatflow/qe/scripts/qe-matdyn.py",
        #"pymatflow/qe/scripts/qe-dynmat.py",
        #"pymatflow/qe/scripts/qe-plotband-for-matdyn.py",
        #"pymatflow/qe/scripts/qe-pp.py",
        #"pymatflow/qe/scripts/qe-phonopy.py",
        #"pymatflow/qe/scripts/qe-molecularpdos.py",
        #"pymatflow/qe/scripts/qe-pw-dielectric.py",
        #"pymatflow/qe/scripts/qe-relax-tetragonal.py",
        #"pymatflow/qe/scripts/qe-relax-cubic.py",
        #"pymatflow/qe/scripts/qe-relax-hexagonal.py",
        #"pymatflow/qe/scripts/qe-get-matdyn-qpoints-from-bands-calc.py",
        #"pymatflow/qe/scripts/qe-get-matdyn-qpoints-from-kpath.py",
        #"pymatflow/qe/post/scripts/post-qe-dos.py",
        #"pymatflow/qe/post/scripts/post-qe-pdos.py",
        #"pymatflow/qe/post/scripts/post-qe-neb.py",
        #"pymatflow/qe/post/scripts/post-qe-scf.py",
        #"pymatflow/qe/post/scripts/post-qe-converge.py",
        #"pymatflow/qe/post/scripts/post-qe-bands.py",
        #"pymatflow/qe/post/scripts/post-qe-matdyn.py",
        #"pymatflow/qe/post/scripts/post-qe-phonopy.py",
        #"pymatflow/siesta/scripts/band_new_scf_siesta.py",
        #"pymatflow/siesta/scripts/charge-density_potentials_siesta.py",
        #"pymatflow/siesta/scripts/chem_analysis_siesta.py",
        #"pymatflow/siesta/scripts/converge_ecut_siesta.py",
        #"pymatflow/siesta/scripts/coord_from_xyz_to_siesta.py",
        #"pymatflow/siesta/scripts/dos-pdos_new_scf_siesta.py",
        #"pymatflow/siesta/scripts/geo_opt_siesta.py",
        #"pymatflow/siesta/scripts/siesta-opt.py",
        #"pymatflow/siesta/scripts/ldos_new_scf_siesta.py",
        #"pymatflow/siesta/scripts/macro_polarization_siesta.py",
        #"pymatflow/siesta/scripts/md_siesta.py",
        #"pymatflow/siesta/scripts/siesta-md.py",
        #"pymatflow/siesta/scripts/net-charge_dipole_elec-field_siesta.py",
        #"pymatflow/siesta/scripts/optical_siesta.py",
        #"pymatflow/siesta/scripts/siesta-phonopy.py",
        #"pymatflow/siesta/scripts/wannier90_siesta.py",
        "pymatflow/siesta/scripts/xyz2fdf.py",
        #"pymatflow/siesta/scripts/siesta-converge-cutoff.py",
        #"pymatflow/siesta/scripts/siesta-scf.py",
        #"pymatflow/siesta/scripts/siesta-scf-restart.py",
        #"pymatflow/siesta/scripts/siesta-ts.py",
        #"pymatflow/siesta/scripts/siesta-phonon.py",
        #"pymatflow/siesta/scripts/siesta-opt-cubic.py",
        #"pymatflow/siesta/scripts/siesta-opt-tetragonal.py",
        #"pymatflow/siesta/scripts/siesta-opt-hexagonal.py",
        #"pymatflow/siesta/post/scripts/post-siesta-pdos.py",
        #"pymatflow/siesta/post/scripts/post-siesta-bands.py",
        #"pymatflow/siesta/post/scripts/post-siesta-converge.py",
        #"pymatflow/siesta/post/scripts/post-siesta-opt.py",
        #"pymatflow/siesta/post/scripts/post-siesta-md.py",
        #"pymatflow/siesta/post/scripts/post-siesta-phonopy.py",
        "pymatflow/scripts/xyz-build-supercell.py",
        "pymatflow/scripts/cluster_sphere.py",
        "pymatflow/scripts/xyz-fix-atoms.py",
        "pymatflow/scripts/seekpath-flow.py",
        "pymatflow/scripts/cif-to-pdb.py",
        "pymatflow/scripts/pdb-to-cif.py",
        "pymatflow/scripts/xyz-modified-to-crystal.py",
        "pymatflow/scripts/kpath-seekpath.py",
        "pymatflow/scripts/kpath-cms-2010.py",
        "pymatflow/scripts/kpoints-mp-suggest.py",
        "pymatflow/scripts/pot-from-xyz-modified.py",
        "pymatflow/scripts/vasp-potcar-from-xyz.py",
        "pymatflow/scripts/xyz-modified-to-poscar.py",
        "pymatflow/scripts/llhpc-env.py",
        "pymatflow/scripts/nebmake.py",
        "pymatflow/scripts/contour-xyz.py",
        "pymatflow/scripts/contour-xyz-diff.py",
        "pymatflow/scripts/ldos-xyz.py",  
        "pymatflow/scripts/xsd-color-bond.py",
        "pymatflow/scripts/stm-vasp.py",
        "pymatflow/scripts/chg-vasp.py",
        "pymatflow/scripts/cube-handle.py",
        "pymatflow/scripts/chg-diff-1d.py",
        "pymatflow/scripts/chgflow.py",
        "pymatflow/scripts/cube-diff-1d.py",
        "pymatflow/scripts/cube-1d.py",
        "pymatflow/scripts/veusz-flow.py",
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
        "pymatflow/flow/scripts/pes-to-img.py",
        #"pymatflow/vasp/scripts/vasp-bands.py",
        #"pymatflow/vasp/scripts/vasp-phono3py.py",
        #"pymatflow/vasp/scripts/vasp-md.py",
        #"pymatflow/vasp/scripts/vasp-phonopy.py",
        #"pymatflow/vasp/scripts/vasp-scf.py",
        #"pymatflow/vasp/scripts/vasp-nscf.py",
        #"pymatflow/vasp/scripts/vasp-converge-encut.py",
        #"pymatflow/vasp/scripts/vasp-converge-sigma.py",
        #"pymatflow/vasp/scripts/vasp-converge-kpoints.py",
        #"pymatflow/vasp/scripts/vasp-opt.py",
        #"pymatflow/vasp/scripts/vasp-neb.py",
        #"pymatflow/vasp/scripts/vasp-dfpt.py",
        #"pymatflow/vasp/scripts/vasp-opt-hexagonal.py",
        #"pymatflow/vasp/scripts/vasp-opt-cubic.py",
        #"pymatflow/vasp/scripts/vasp-opt-tetragonal.py",
        #"pymatflow/vasp/scripts/vasp-efficiency-test.py",
        #"pymatflow/vasp/scripts/vasp-phonon.py",
        #"pymatflow/vasp/scripts/vasp-berry.py",
        "pymatflow/vasp/scripts/vasp-time.py",
        "pymatflow/vasp/scripts/vasp-grep-energy.py",
        "pymatflow/vasp/scripts/xdatcar-to-xyz.py",
        #"pymatflow/vasp/scripts/split-poscars.py",
        #"pymatflow/vasp/post/scripts/post-vasp-opt.py",
        #"pymatflow/vasp/post/scripts/post-vasp-converge.py",
        #"pymatflow/vasp/post/scripts/post-vasp-scf.py",
        #"pymatflow/vasp/post/scripts/post-vasp-phonopy.py",
        #"pymatflow/vasp/post/scripts/post-vasp-bands-p4vasp.py",
        #"pymatflow/vasp/post/scripts/post-vasp-opt-dev.py",
        #"pymatflow/vasp/post/scripts/post-vasp-scf-dev.py",
        #"pymatflow/vasp/post/scripts/post-vasp-phonon.py",
        #"pymatflow/vasp/post/scripts/post-vasp-bands.py",
        #"pymatflow/vasp/post/scripts/post-vasp-pdos.py",
        #"pymatflow/vasp/post/scripts/post-vasp-neb-vtst.py",
        #"pymatflow/vasp/post/scripts/procar-vasp.py",
        "pymatflow/vasp/post/scripts/vasp-tensor-piezo.py",
        "pymatflow/octopus/scripts/octflow.py",
        ],
    entry_points = {
        'console_scripts': [
            'matflow = pymatflow.cmd.matflow:main',
            'mflow = pymatflow.cmd.matflow:main',
            'postflow = pymatflow.cmd.postflow:main',
            'pflow = pymatflow.cmd.postflow:main',
            'structflow = pymatflow.cmd.structflow:main',
            'sflow = pymatflow.cmd.structflow:main',
            'lmpflow = pymatflow.cmd.lmpflow:main',
            'lflow = pymatflow.cmd.lmpflow:main',
            "psoflow = pymatflow.flow.calypso.psoflow:main",
            "octoflow = pymatflow.octopus.octoflow:main",
            "excitingflow = pymatflow.exciting.excitingflow:main",
            "elkflow = pymatflow.elk.elkflow:main",
            "xtbflow = pymatflow.xtb.xtbflow:main",
        ]
    },
)

