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



def pybind11_build_ext(builder, ext):
    extdir = os.path.abspath(os.path.dirname(builder.get_ext_fullpath(ext.name)))
    #extdir = os.path.join(extdir, "pymatflow/cpp")
    # required for auto-detection & inclusion of auxiliary "native" libs
    if not extdir.endswith(os.path.sep):
        extdir += os.path.sep
    debug = int(os.environ.get("DEBUG", 0)) if builder.debug is None else builder.debug
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
        "-DEXAMPLE_VERSION_INFO={}".format(builder.distribution.get_version())
    ]
    if builder.compiler.compiler_type != "msvc":
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
            cmake_args += ["-A", PLAT_TO_CMAKE[builder.plat_name]]
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
        if hasattr(builder, "parallel") and builder.parallel:
            # CMake 3.12+ only.
            build_args += ["-j{}".format(builder.parallel)]
    if not os.path.exists(builder.build_temp):
        os.makedirs(builder.build_temp)
    subprocess.check_call(
        ["cmake", ext.sourcedir] + cmake_args, cwd=builder.build_temp
    )
    subprocess.check_call(
        ["cmake", "--build", "."] + build_args, cwd=builder.build_temp
    )
    #
    #if ext.name in ["pyaskit"]:
    #    os.system("mv _skbuild/linux-x86_64-*/setuptools/lib.linux-x86_64-*/%s.*.so _skbuild/linux-x86_64-*/setuptools/lib.linux-x86_64-*/pymatflow/cpp/" % ext.name)

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

# --------------------------------
# fortran support
# --------------------------------
from Cython.Build import cythonize
ext_modules_fortran = cythonize(
    [
        Extension(
            'cube_handle',
            sources=[
                "fortran/pyx/cube_handle.pyx"
            ],
            # other compile args for gcc
            extra_compile_args=['-fPIC', '-O3'],
            # other files to link to
            extra_link_args=[
                #'fortran/c_binding/libaskitf-c-binding.a',
                'fortran/build/libaskitf-c-binding.a',
                #'fortran/src/libaskitf.a',
                'fortran/build/libaskitf.a',
                "-lgfortran",
                "-fopenmp" # important for OpenMP dependency
            ],
            language="c"
        ),
    ],
    language_level="3"
)

ext_modules_all = []
ext_modules_all += ext_modules_fortran

ext_modules_fortran_names = []
for ext in ext_modules_fortran:
    ext_modules_fortran_names.append(ext.name)

def fortran_build_ext(builder, ext):
    builder.build_extension(ext)
    # TODO: more strong code here
    os.system("mv _skbuild/linux-x86_64-*/setuptools/lib.linux-x86_64-*/%s.*.so _skbuild/linux-x86_64-*/setuptools/lib.linux-x86_64-*/pymatflow/fortran/" % ext.name)


class CustomBuildExt(build_ext):
    """
    """
    def run(self):
        #os.system("make -C ./fortran/src")
        #os.system("make -C ./fortran/cmd")
        #os.system("make -C ./fortran/c_binding")
        os.system("mkdir -p ./fortran/build")
        os.chdir("./fortran/build")
        os.system("cmake ..; make")
        os.chdir("../../")
        super().run()

    def build_extension(self, ext):
        if ext.name in ext_modules_fortran_names:
            fortran_build_ext(super(), ext)
        else:
            pybind11_build_ext(self, ext)


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
    cmdclass = {
        'clean': CleanCommand,
        'build_ext': CustomBuildExt
    },    
    ext_modules=[
        # --------------------------------------------
        # cpp extension using CustomBuildExt (working)
        # --------------------------------------------
        #CMakeExtension(name="pyaskit", sourcedir="cpp")
        # --------------------------------------------
    ]+ext_modules_fortran,
    # ------------------------------------------
    # cpp extension using scikit-build (working)
    # ------------------------------------------
    cmake_source_dir="cpp", # where CMakeLists.txt exists
    cmake_install_dir="pymatflow/cpp", # from pymatflow.cpp import pyaskit
    cmake_args=["-DCMAKE_BUILD_TYPE=Debug"],
    # ------------------------------------------
    scripts = [
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
        "pymatflow/vasp/scripts/vasp-time.py",
        "pymatflow/vasp/scripts/vasp-grep-energy.py",
        "pymatflow/vasp/scripts/xdatcar-to-xyz.py",
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

