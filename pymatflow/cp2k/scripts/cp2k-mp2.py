#!/usr/bin/env python
# _*_ coding: utf-8 _*_

import os
import argparse

from pymatflow.cp2k.mp2 import static_mp2_run
from pymatflow.remote.server import server_handle
"""
usage:
"""


params = {}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--directory", type=str, default="tmp-cp2k-static-mp2",
            help="directory where the calculation happens")

    parser.add_argument("-f", "--file", type=str,
            help="the xyz file containing the structure to be simulated")


    parser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    parser.add_argument("--mpi", type=str, default="",
            help="mpi command: like --mpi='mpirun -np 4'")

    # ------------------------------------------------------------------
    #                    force_eval/dft related parameters
    # ------------------------------------------------------------------

    parser.add_argument("--qs-method", type=str, default="gpw",
            choices=["am1", "dftb", "gapw", "gapw_xc", "gpw", "lrigpw", "mndo", "mndod",
                "ofgpw", "pdg", "pm3", "pm6", "pm6-fm", "pnnl", "rigpw", "rm1"],
            help="dft-qs-method: specify the electronic structure method")

    parser.add_argument("--eps-scf", type=float, default=1.0e-6,
            help="dft-scf-eps_scf")

    parser.add_argument("--xc-functional", type=str, default="pbe",
            help="dft-xc-xc_functional: LYP, PADE, PBE, PW92, TPSS, XGGA, XWPBE, etc.")

    parser.add_argument("--cutoff", type=int, default=100,
            help="CUTOFF, default value: 100 Ry")

    parser.add_argument("--rel-cutoff", type=int, default=60,
            help="REL_CUTOFF, default value: 60 Ry")

    parser.add_argument("-k", "--kpoints-scheme", type=str,
            default="GAMMA",
            help="DFT-KPOINTS-SCHEME(str): can be NONE, GAMMA, MONKHORST-PACK, MACDONALD, GENERAL. when you set MONKHORST-PACK, you should also add the three integers like 'monkhorst-pack 3 3 3'")

    parser.add_argument("--diag", type=str, default="TRUE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing tranditional diagonalization for SCF")

    parser.add_argument("--ot", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing orbital transformation for SCF")

    parser.add_argument("--alpha", type=float, default=0.4,
            help="DFT-SCF-MIXING-ALPHA")

    parser.add_argument("--smear", type=str, default="FALSE",
            #choices=["TRUE", "FALSE", "true", "false"],
            help="switch on or off smearing for occupation")

    parser.add_argument("--smear-method", type=str, default="FERMI_DIRAC",
            help="smearing type: FERMI_DIRAC, ENERGY_WINDOW")

    parser.add_argument("--added-mos", type=int, default=0,
            help="ADDED_MOS for SCF")


    parser.add_argument("--electronic-temp", type=float, default=300,
            help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR")

    parser.add_argument("--window-size", type=float, default=0,
            help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing")

    parser.add_argument("--ls-scf", type=str, default="false",
            #choices=["true", "false", "true", "false"],
            help="dft-ls_scf: use linear scaling scf method")

    # vdw correction related
    parser.add_argument("--vdw", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether to use VDW correction")

    parser.add_argument("--vdw-potential-type", type=str, default="PAIR_POTENTIAL",
            choices=["PAIR_POTENTIAL", "NON_LOCAL", "NONE"],
            help="DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE: PAIR_POTENTIAL, NON_LOCAL")

    parser.add_argument("--pair-type", type=str, default="DFTD3",
            choices=["DFTD2", "DFTD3", "DFTD3(BJ)"],
            help="VDW PAIR_POTENTIAL type: DFTD2, DFTD3, DFTD3(BJ)")

    parser.add_argument("--r-cutoff", type=float, default=1.05835442E+001,
            help="DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL: Range of potential. The cutoff will be 2 times this value")

    parser.add_argument("-p", "--printout-option", nargs="+", type=int,
            default=[],
            choices = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
            help=
            """
            Properties printout option, you can also activate multiple prinout-option at the same time.
            1: printout pdos
            2: printout band
            3: printout electron densities
            4: printout electron local function(ELF)
            5: printout molecular orbitals
            6: printout molecular orbital cube files
            7: printout mulliken populaltion analysis
            8: printout cubes for generation of STM images
            9: printout cube file with total density(electrons+atomic core)
           10: printout v_hartree_cube
           11: printout v_xc_cube
           12: printout xray_diffraction_spectrum
           13: request a RESP fit of charges.
           default is no printout of these properties.
           """)

    parser.add_argument("--dft-print-elf-cube-stride", type=int, nargs="+",
            default=[1, 1, 1],
            help="DFT-PRINT-ELF_CUBE-STRIDE")

    parser.add_argument("--dft-print-e-density-cube-stride", type=int, nargs="+",
            default=[1, 1, 1],
            help="DFT-PRINT-E_DENSITY_CUBE-STRIDE")

    # ------------------------------------------------------------------
    #                    force_eval/properties related parameters
    # ------------------------------------------------------------------

    parser.add_argument("--properties-resp-slab-sampling-range", type=float, nargs="+",
            default=[0.3, 3.0],
            help="PROPERTIES-RESP-SLAB_SAMPLING-RANGE.")

    parser.add_argument("--properties-resp-slab-sampling-surf-direction", type=str, default="Z",
            choices=["X", "Y", "Z", "x", "y", "z", "-X", "-Y", "-Z", "-x", "-y", "-z"],
            help="PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION.")

    parser.add_argument("--properties-resp-slab-sampling-atom-list", type=int, nargs="+",
            default=[1],
            help="PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST")


    # -----------------------------------------------------------
    #            ATOM related parameters
    # -----------------------------------------------------------


    # -----------------------------------------------------------------
    #                      for server handling
    # -----------------------------------------------------------------
    parser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")
    parser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")
    parser.add_argument("--jobname", type=str, default="geo-opt",
            help="jobname on the pbs server")
    parser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")
    parser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # ==========================================================
    # transfer parameters from the arg parser to opt_run setting
    # ==========================================================
    args = parser.parse_args()

    params["FORCE_EVAL-DFT-LS_SCF"] = args.ls_scf
    params["FORCE_EVAL-DFT-QS-METHOD"] = args.qs_method
    params["FORCE_EVAL-DFT-MGRID-CUTOFF"] = args.cutoff
    params["FORCE_EVAL-DFT-MGRID-REL_CUTOFF"] = args.rel_cutoff
    params["FORCE_EVAL-DFT-XC-XC_FUNCTIONAL"] = args.xc_functional
    params["FORCE_EVAL-DFT-SCF-EPS_SCF"] = args.eps_scf
    params["FORCE_EVAL-DFT-SCF-ADDED_MOS"] = args.added_mos
    params["FORCE_EVAL-DFT-SCF-SMEAR"] = args.smear
    params["FORCE_EVAL-DFT-SCF-SMEAR-METHOD"] = args.smear_method
    params["FORCE_EVAL-DFT-SCF-SMEAR-ELECTRONIC_TEMPERATURE"] = args.electronic_temp
    params["FORCE_EVAL-DFT-SCF-SMEAR-WINDOW_SIZE"] = args.window_size
    params["FORCE_EVAL-DFT-SCF-DIAGONALIZATION"] = args.diag
    params["FORCE_EVAL-DFT-SCF-OT"] = args.ot
    params["FORCE_EVAL-DFT-SCF-MIXING-ALPHA"] = args.alpha
    params["FORCE_EVAL-DFT-KPOINTS-SCHEME"] = args.kpoints_scheme

    params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL"] = args.vdw
    params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE"] = args.vdw_potential_type
    parmas["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-TYPE"] = args.pair_type
    params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR-POTENTIAL-R_CUTOFF"] = args.r_cutoff


    params["FORCE_EVAL-DFT-PRINT-ELF_CUBE-STRIDE"] = args.dft_print_elf_cube_stride
    params["FORCE_EVAL-DFT-PRINT-E_DENSITY_CUBE-STRIDE"] = args.dft_print_e_density_cube_stride


    params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-RANGE"] = args.properties_resp_slab_sampling_range
    params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION"] = args.properties_resp_slab_sampling_surf_direction
    params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST"] = args.properties_resp_slab_sampling_atom_list


    task = static_mp2_run()
    task.get_xyz(args.file)
    task.set_params(params=params)
    task.set_printout(args.printout_option)
    task.scf_mp2(directory=args.directory, mpi=args.mpi, runopt=args.runopt, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)

    server_handle(auto=args.auto, directory=args.directory, jobfilebase="static-scf-mp2", server=args.server)
