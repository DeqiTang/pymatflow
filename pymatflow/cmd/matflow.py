#!/usr/bin/env python

import os
import sys
import argparse



def get_kpath(kpath_manual=None, kpath_file=None):
    """
    :param kpath_manual: manual input kpath like --kpath '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '
    :param kpath_file: manual input kpath read from the file
    :return kpath or None(when kpath_manual and kpath_file are both None)
    """
    # dealing with standard kpath
    kpath = None
    if kpath_manual != None:
        # kpath from script argument args.kpath
        kpath = []
        for kpoint in kpath_manual:
            if kpoint.split()[4] != "|":
                kpath.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    int(kpoint.split()[4]),
                    ])
            elif kpoint.split()[4] == "|":
                kpath.append([
                    float(kpoint.split()[0]),
                    float(kpoint.split()[1]),
                    float(kpoint.split()[2]),
                    kpoint.split()[3].upper(),
                    "|",
                    ])
        return kpath
    elif kpath_file != None:
        # kpath read from file specified by kpath_file
        # file is in format like this
        """
        5
        0.0 0.0 0.0 #GAMMA 15
        x.x x.x x.x #XXX |
        x.x x.x x.x #XXX 10
        x.x x.x x.x #XXX 15
        x.x x.x x.x #XXX 20
        """
        # if there is a '|' behind the label it means the path is
        # broken after that point!!!
        kpath = []
        with open(kpath_file, 'r') as fin:
            lines = fin.readlines()
        nk = int(lines[0])
        for i in range(nk):
            if lines[i+1].split("\n")[0].split()[4] != "|":
                kpath.append([
                    float(lines[i+1].split()[0]),
                    float(lines[i+1].split()[1]),
                    float(lines[i+1].split()[2]),
                    lines[i+1].split()[3].split("#")[1].upper(),
                    int(lines[i+1].split()[4]),
                    ])
            elif lines[i+1].split("\n")[0].split()[4] == "|":
                kpath.append([
                    float(lines[i+1].split()[0]),
                    float(lines[i+1].split()[1]),
                    float(lines[i+1].split()[2]),
                    lines[i+1].split()[3].split("#")[1].upper(),
                    '|',
                    ])
        return kpath
    else:
        pass
    # -------------------------------------------------------------------



def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="driver", title="subcommands", description="choose one and only one calculator")

    # --------------------------------------------------------------------------
    # Abinit
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("abinit", help="using abinit as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default="static",
            choices=[0, 1, 2, 3],
            help="choices of runtype. 0->static_run; 1->optimization; 2->dfpt-elastic-piezo-dielec")

    subparser.add_argument("--chkprim", type=int, default=1,
            choices=[0, 1],
            help="check whether the input cell is primitive. if your cell is not primitive, set chkprim to 0. for more information, refer to https://docs.abinit.org/variables/gstate/#chkprim")


    subparser.add_argument("--properties", nargs="+", type=int,
            default=[],
            help="options for properties calculation")

    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    subparser.add_argument("--kpath-file", type=str,
            help="file to read the kpath for band structure calculation")

    subparser.add_argument("--iscf", type=int, default=7,
            choices=[0, 1, 2, 3, 4, 5, 7, 12, 13, 14, 15, 17],
            help="set scf or nscf type. for more information, refer to https://docs.abinit.org/variables/basic/#iscf")

    subparser.add_argument("--ecut", type=int, default=15,
            help="Kinetic energy cutoff for wave functions in unit of Hartree, default value: 15 Hartree. for more information, refer to https://docs.abinit.org/variables/basic/#ecut")

    subparser.add_argument("--ixc", type=int, default=11,
            choices=[1, 2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 40, 41, 42],
            help="type of exchage-correlation functional. for more information, refer to https://docs.abinit.org/variables/basic/#ixc")

    subparser.add_argument("--kptopt", type=int, default=1,
            choices=[1],
            help="Kpoints Generation scheme option: 0, 1, 2, 3, 4 or a negative value. for more information, refer to https://docs.abinit.org/variables/basic/#kptopt")

    subparser.add_argument("--ngkpt", nargs="+", type=int,
            default=[1, 1, 1],
            help="number of grid points for kpoints generation. for more information, refer to https://docs.abinit.org/variables/basic/#ngkpt")

    # electron occupation
    subparser.add_argument("--occopt", type=int, default=3,
            choices=[0, 1, 2, 3, 4, 5, 6, 7],
            help="Controls how input parameters nband, occ, and wtk are handled. for more information, refer to https://docs.abinit.org/variables/basic/#occopt")

    subparser.add_argument("--nband", type=int, nargs="+", default=None,
            help="Gives number of bands, occupied plus possibly unoccupied, for which wavefunctions are being computed along with eigenvalues. for more information, refer to https://docs.abinit.org/variables/basic/#nband")

    subparser.add_argument("--occ", nargs="+", type=float, default=None,
            help="Gives occupation numbers for all bands in the problem. Needed if occopt == 0 or occopt == 2. Ignored otherwise. Also ignored when iscf = -2. refer to https://docs.abinit.org/variables/gstate/#occ")

    # magnetic related parameters
    subparser.add_argument("--nsppol", type=int, default=None,
            choices=[1, 2],
            help="Give the number of INDEPENDENT spin polarisations, for which there are non- related wavefunctions. Can take the values 1 or 2. for more information, refer to https://docs.abinit.org/variables/basic/#nsppol")

    # vdw related parameters
    subparser.add_argument("--vdw-xc", type=int, default=None,
            choices=[0, 1, 2, 5, 6, 7, 10, 11, 14],
            help="Van Der Waals exchange-correlation functional. 0: no correction, 1: vdW-DF1, 2: vdW-DF2, 5: DFT-D2, 6: DFT-D3, 7: DFT-D3(BJ). for more information, refer to https://docs.abinit.org/variables/vdw/#vdw_xc")

    subparser.add_argument("--vdw-tol", type=float,
            default=None,
            help="Van Der Waals tolerance, only work when vdw_xc == 5 or 6 or 7. to be included in the potential a pair of atom must have contribution to the energy larger than vdw_tol. default value is 1.0e-10. fore more information, refer to https://docs.abinit.org/variables/vdw/#vdw_tol")

    subparser.add_argument("--prtden", type=int ,default=1,
            choices=[0, 1],
            help="print the density. for more information, refer to https://docs.abinit.org/variables/files/#prtden")

    subparser.add_argument("--prtdos", type=int, default=None,
            choices=[0, 1, 2, 3],
            help="can be 0, 1, 2, 3. for more information, refer to https://docs.abinit.org/variables/files/#prtdos")


    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory to do the calculation")

    # structure file: either xyz or cif. they are exclusive
    # actually this can be put in the main subparser, but it will make the command not like git sub-cmmand
    # so we put them in every subsubparser
    structfile = subparser.add_mutually_exclusive_group(required=True) # at leaset one of cif and xyz is provided
    # argparse will make sure only one of argument in structfile(xyz, cif) appear on command line
    structfile.add_argument("--xyz", type=str, default=None,
            help="The xyz structure file with second line specifying cell parameters")

    structfile.add_argument("--cif", type=str, default=None,
            help="The cif structure file")

    # run option
    subparser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    subparser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")


    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    subparser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    subparser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    subparser.add_argument("--jobname", type=str, default="matflow-job",
            help="jobname on the pbs server")

    subparser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    subparser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # --------------------------------------------------------------------------
    # CP2K
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("cp2k", help="using cp2k as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default="static",
            choices=[0, 1, 2, 3],
            help="choices of runtype. 0->static_run; 1->optimization;")

    # force_eval/dft related parameters

    subparser.add_argument("--qs-method", type=str, default="gpw",
            choices=["am1", "dftb", "gapw", "gapw_xc", "gpw", "lrigpw", "mndo", "mndod",
                "ofgpw", "pdg", "pm3", "pm6", "pm6-fm", "pnnl", "rigpw", "rm1"],
            help="dft-qs-method: specify the electronic structure method")

    subparser.add_argument("--eps-scf", type=float, default=1.0e-6,
            help="dft-scf-eps_scf")

    subparser.add_argument("--xc-functional", type=str, default="pbe",
            help="dft-xc-xc_functional: LYP, PADE, PBE, PW92, TPSS, XGGA, XWPBE, etc.")

    subparser.add_argument("--cutoff", type=int, default=100,
            help="CUTOFF, default value: 100 Ry")

    subparser.add_argument("--rel-cutoff", type=int, default=60,
            help="REL_CUTOFF, default value: 60 Ry")

    subparser.add_argument("-k", "--kpoints-scheme", type=str,
            default="GAMMA",
            help="DFT-KPOINTS-SCHEME(str): can be NONE, GAMMA, MONKHORST-PACK, MACDONALD, GENERAL. when you set MONKHORST-PACK, you should also add the three integers like 'monkhorst-pack 3 3 3'")

    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath for band structure calculation")

    subparser.add_argument("--kpath-file", type=str,
            help="file to read the kpath for band structure calculation")

    subparser.add_argument("--diag", type=str, default="TRUE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing tranditional diagonalization for SCF")

    subparser.add_argument("--ot", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether choosing orbital transformation for SCF")

    subparser.add_argument("--alpha", type=float, default=0.4,
            help="DFT-SCF-MIXING-ALPHA")

    subparser.add_argument("--smear", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="switch on or off smearing for occupation")

    subparser.add_argument("--smear-method", type=str, default="FERMI_DIRAC",
            help="smearing type: FERMI_DIRAC, ENERGY_WINDOW")

    subparser.add_argument("--added-mos", type=int, default=0,
            help="ADDED_MOS for SCF")

    subparser.add_argument("--electronic-temp", type=float, default=300,
            help="ELECTRON_TEMPERATURE for FERMI_DIRAC SMEAR")

    subparser.add_argument("--window-size", type=float, default=0,
            help="Size of the energy window centred at the Fermi level for ENERGY_WINDOW type smearing")

    subparser.add_argument("--ls-scf", type=str, default="false",
            choices=["true", "false", "true", "false"],
            help="dft-ls_scf: use linear scaling scf method")

    # vdw correction related
    subparser.add_argument("--usevdw", type=str, default="FALSE",
            choices=["TRUE", "FALSE", "true", "false"],
            help="whether to use VDW correction")

    subparser.add_argument("--vdw-potential-type", type=str, default="PAIR_POTENTIAL",
            choices=["PAIR_POTENTIAL", "NON_LOCAL", "NONE"],
            help="DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE: PAIR_POTENTIAL, NON_LOCAL")

    subparser.add_argument("--pair-type", type=str, default="DFTD3",
            choices=["DFTD2", "DFTD3", "DFTD3(BJ)"],
            help="VDW PAIR_POTENTIAL type: DFTD2, DFTD3, DFTD3(BJ)")

    subparser.add_argument("--r-cutoff", type=float, default=1.05835442E+001,
            help="DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL: Range of potential. The cutoff will be 2 times this value")

    subparser.add_argument("-p", "--printout-option", nargs="+", type=int,
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

    subparser.add_argument("--dft-print-elf-cube-stride", type=int, nargs="+",
            default=[1, 1, 1],
            help="DFT-PRINT-ELF_CUBE-STRIDE")

    subparser.add_argument("--dft-print-e-density-cube-stride", type=int, nargs="+",
            default=[1, 1, 1],
            help="DFT-PRINT-E_DENSITY_CUBE-STRIDE")

    # ------------------------------------------------------------------
    #                    force_eval/properties related parameters
    # ------------------------------------------------------------------

    subparser.add_argument("--properties-resp-slab-sampling-range", type=float, nargs="+",
            default=[0.3, 3.0],
            help="PROPERTIES-RESP-SLAB_SAMPLING-RANGE.")

    subparser.add_argument("--properties-resp-slab-sampling-surf-direction", type=str, default="Z",
            choices=["X", "Y", "Z", "x", "y", "z", "-X", "-Y", "-Z", "-x", "-y", "-z"],
            help="PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION.")

    subparser.add_argument("--properties-resp-slab-sampling-atom-list", type=int, nargs="+",
            default=[1],
            help="PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST")

    # structure file: either xyz or cif. they are exclusive
    # actually this can be put in the main subparser, but it will make the command not like git sub-cmmand
    # so we put them in every subsubparser
    structfile = subparser.add_mutually_exclusive_group(required=True) # at leaset one of cif and xyz is provided
    # argparse will make sure only one of argument in structfile(xyz, cif) appear on command line
    structfile.add_argument("--xyz", type=str, default=None,
            help="The xyz structure file with second line specifying cell parameters")

    structfile.add_argument("--cif", type=str, default=None,
            help="The cif structure file")

    # run option
    subparser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    subparser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")


    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    subparser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    subparser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    subparser.add_argument("--jobname", type=str, default="matflow-job",
            help="jobname on the pbs server")

    subparser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    subparser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")


    # --------------------------------------------------------------------------
    # Quantum ESPRESSO
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("qe", help="using quantum espresso as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1, 2, 3, 4],
            help="choices of runtype. 0->static_run; 1->optimization;")

    subparser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")


    # -------------------------------------------------------------------
    #                       scf related parameters
    # -------------------------------------------------------------------
    subparser.add_argument("--ecutwfc",
            type=int, default=100)

    subparser.add_argument("--ecutrho", type=int, default=None,
            help="Kinetic energy cutoff for charge density and potential in unit of Rydberg, default value: None")

    subparser.add_argument("--kpoints-option", type=str, default="automatic",
            choices=["automatic", "gamma", "crystal_b"],
            help="Kpoints generation scheme option for the SCF or non-SCF calculation")

    subparser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[1, 1, 1, 0, 0, 0],
            help="Monkhorst-Pack kpoint grid, in format like --kpoints-mp 1 1 1 0 0 0")

    subparser.add_argument("--kpath-manual", type=str, nargs="+", default=None,
            help="manual input kpath in crystal_b, like --kpath-manual '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '")

    subparser.add_argument("--kpath-file", type=str,
            help="manual input kpath in crystal_b read from the file")


    subparser.add_argument("--conv-thr", type=float, default=1.0e-6,
            help="the conv_thr for scf, when doing geometric optimization better use a strict covnergec for scf")

    subparser.add_argument("--occupations", type=str, default="smearing",
            choices=["smearing", "tetrahedra", "tetrahedra_lin", "tetrahedra_opt", "fixed", "from_input"],
            help="Occupation method for the calculation.")

    subparser.add_argument("--smearing", type=str, default="gaussian",
            choices=["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"],
            help="Smearing type for occupations by smearing, default is gaussian in this script")

    subparser.add_argument("--degauss", type=float, default=0.001,
            help="Value of the gaussian spreading (Ry) for brillouin-zone integration in metals.(defualt: 0.001 Ry)")

    subparser.add_argument("--nbnd", type=int, default=None,
            help="Number of electronic states (bands) to be calculated")

    subparser.add_argument("--tstress", type=str, default=".false.",
            choices=[".true.", ".false."],
            help="calculate stress. default=.false.")

    subparser.add_argument("--vdw-corr", help="vdw_corr = dft-d, dft-d3, ts-vdw, xdm", type=str, default="none")

    # magnetic related parameters
    subparser.add_argument("--nspin", type=int, default=None,
            choices=[1, 2],
            help="choose either 1 or 2, and 4 should not be used as suggested by pwscf official documentation.")

    subparser.add_argument("--starting-magnetization", type=float, nargs="+", default=None,
            help="starting_magnetization(i), i=1,ntyp -> Starting spin polarization on atomic type i in a spin polarized calculation. Values range between -1 (all spins down for the valence electrons of atom type i) to 1 (all spins up).")

    subparser.add_argument("--noncolin", type=str, default=None,
            choices=[".true.", ".false."],
            help="if .true. the program will perform a noncollinear calculation.")

    # ATOMIC_FORCES
    subparser.add_argument("--pressure", type=float, default=None,
            help="specify pressure acting on system in unit of Pa")
    subparser.add_argument("--pressuredir", type=str, default=None,
            choices=["x", "y", "z"],
            help="specify direction of pressure acting on system.")

    # -------------------------------------------------------------------
    #               geometric optimization related parameters
    # -------------------------------------------------------------------
    subparser.add_argument("--etot-conv-thr",
            type=float, default=1.0e-4,
            help="convergence threshold of energy for geometric optimization")

    subparser.add_argument("--forc-conv-thr",
            type=float, default=1.0e-3,
            help="convergence threshold for force in optimization,(usually it is more important than energy)")

    subparser.add_argument("--nstep",
            type=int, default=50,
            help="maximum ion steps for geometric optimization")

    subparser.add_argument("--cell-dofree", type=str, default=None,
            choices=['all', 'ibrav', 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz', 'shape', 'volume', '2Dxy', '2Dshape', 'epitaxial_ab', 'epitaxial_ac', 'epitaxial_bc'],
            help="cell_dofree for &cell/")

    # structure file: either xyz or cif. they are exclusive
    # actually this can be put in the main subparser, but it will make the command not like git sub-cmmand
    # so we put them in every subsubparser
    structfile = subparser.add_mutually_exclusive_group(required=True) # at leaset one of cif and xyz is provided
    # argparse will make sure only one of argument in structfile(xyz, cif) appear on command line
    structfile.add_argument("--xyz", type=str, default=None,
            help="The xyz structure file with second line specifying cell parameters")

    structfile.add_argument("--cif", type=str, default=None,
            help="The cif structure file")

    subparser.add_argument("--images", type=str, nargs="+",
            help="the image xyz file(--images first.xyz imtermediate-1.xyz intermediate-2.xyz ... last.xyz)")

    # run option
    subparser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    subparser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")

    # params for neb namelist &path
    subparser.add_argument("--string-method", type=str, default="neb",
            help="string_method")

    subparser.add_argument("--nstep-path", type=int, default=100,
            help="nstep_path")

    subparser.add_argument("--opt-scheme", type=str, default="broyden",
            help="Specify the type of optimization scheme(sd, broyden, broyden2, quick-min, langevin)")

    subparser.add_argument("--num-of-images", type=int, default=5,
            help="number of total images(including the initial and final image). about how to set proper number of images: usually the inter-image distance between 1~2Bohr is OK")

    subparser.add_argument("--k-max", type=float, default=0.3e0,
            help="Set them to use a Variable Elastic Constants scheme elastic constants are in the range [ k_min, k_max  ], this is useful to rise the resolution around the saddle point")

    subparser.add_argument("--k-min", type=float, default=0.2e0,
            help="Set them to use a Variable Elastic Constants scheme elastic constants are in the range [ k_min, k_max  ], this is useful to rise the resolution around the saddle point")

    subparser.add_argument("--ci-scheme", type=str, default="auto",
            help="Specify the type of Climbing Image scheme(no-CI, auto, manual)")

    subparser.add_argument("--path_thr", type=float, default=0.05,
            help="path_thr")

    subparser.add_argument("--ds", type=float, default=1.e0, help="Optimisation step length ( Hartree atomic units )")

    subparser.add_argument("--first-last-opt", type=bool, default=False,
            help="whether to optimize the first and last image")


    # for phx
    # --------------------------------------------------------------
    parser.add_argument("--tr2-ph", type=float, default=1.0e-14,
            help="threshold for self-consistency.")

    parser.add_argument("--nq", type=int, nargs="+",
            default=[0, 0, 0],
            help="set value of nq1 nq2 nq3.")

    parser.add_argument("--epsil", type=str, default=None,
            choices=[".true.", ".false."],
            help="set epsil in inputph")

    parser.add_argument("--lraman", type=str, default=None,
            choices=["true", "false"],
            help="set lraman, can be 'true' or 'false' only. default is None which means 'false' in real world.")



    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    subparser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    subparser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    subparser.add_argument("--jobname", type=str, default="matflow-job",
            help="jobname on the pbs server")

    subparser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    subparser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")



    # --------------------------------------------------------------------------
    # SIESTA
    # --------------------------------------------------------------------------
    subparser = subparsers.add_parser("siesta", help="using siesta as calculator")

    subparser.add_argument("-r", "--runtype", type=int, default=0,
            choices=[0, 1],
            help="choices of runtype. 0->static_run; 1->optimization;")

    parser.add_argument("-d", "--directory", type=str, default="matflow-running",
            help="Directory for the running.")

    # --------------------------------------------------------------------------
    subparser.add_argument("--meshcutoff", type=int, default=200,
            help="MeshCutoff (Ry)")

    subparser.add_argument("--solution-method", type=str, default="diagon",
            choices=["diagon", "OMM", "OrderN", "PEXSI"],
            help="SolutionMethod(diagon, OMM, OrderN, PEXSI)")

    subparser.add_argument("--functional", type=str, default="GGA",
            help="XC.functional")

    subparser.add_argument("--authors", type=str, default="PBE",
            help="XC.authors")

    subparser.add_argument("--tolerance", type=float, default=1.0e-6,
            help="DM.Tolerance")

    subparser.add_argument("--numberpulay", type=int, default=8,
            help="DM.NumberPulay")

    subparser.add_argument("--mixing", type=float, default=0.1,
            help="DM.MixingWeight")

    subparser.add_argument("--kpoints-mp", type=int, nargs="+",
            default=[3, 3, 3],
            help="set kpoints like '3 3 3'")

    subparser.add_argument("--occupation", type=str, default="FD",
            choices=["FD", "MP"],
            help="OccupationFunction(FD or MP)")

    subparser.add_argument("--electronic-temperature", type=int, default=300,
            help="Electronic Temperature")


    # ------------------------------
    # properties related parameter
    # ------------------------------
    subparser.add_argument("-p", "--properties" ,nargs="+", type=int, default=[],
            help="Option for properties calculation")

    subparser.add_argument("--pdos-block", type=float, nargs="+",
            default=[-20, 10, 0.2, 500])
    #------------------------------------------------------------------------------------------------
    subparser.add_argument("--kpath", type=str, nargs="+", default=None,
            help="manual input kpath in bandlines modes, like --kpath '0.000000 0.000000 0.000000 GAMMA 5' '0.500000 0.000000 0.000000 X 5' '0.0000 0.000 0.50000 A |' '0.5 0.5 0.5 R '")

    subparser.add_argument("--kpath-manual-file", type=str,
            help="manual input kpath in bandline mode  read from the file")
    #------------------------------------------------------------------------------------------------
    subparser.add_argument("--polarization-grids", nargs="+", type=str,
            default=["10 3 3 no", "2 20 2 no", "4 4 15 no"],
            help="PolarizationGrids")

    subparser.add_argument("--external-electric-field", nargs="+", type=float,
            default=[0.0, 0.0, 0.5],
            help="External Electric field")

    subparser.add_argument("--optical-energy-minimum", type=float,
            default=0.0,
            help="Optical.Energy.Minimum")

    subparser.add_argument("--optical-energy-maximum", type=float,
            default=10.0,
            help="Optical.Energy.Maximum")

    subparser.add_argument("--optical-broaden", type=float,
            default=0.0,
            help="Optical.Broaden")

    subparser.add_argument("--optical-scissor", type=float,
            default=0.0,
            help="Optical.Scissor")

    subparser.add_argument("--optical-mesh", nargs="+", type=int,
            default=[5, 5, 5],
            help="Optical.Mesh")

    subparser.add_argument("--optical-polarization-type", type=str,
            default="unpolarized",
            help="Optical.PolarizationType")

    subparser.add_argument("--optical-vector", nargs="+", type=float,
            default=[1.0, 0.0, 0.5],
            help="Optical.Vector")

    subparser.add_argument("--wannier90-unkgrid", nargs="+", type=int,
            default=[10, 10, 10],
            help="Siesta2Wannier90.UnkGrid[1-3]")

    # structure file: either xyz or cif. they are exclusive
    # actually this can be put in the main subparser, but it will make the command not like git sub-cmmand
    # so we put them in every subsubparser
    structfile = subparser.add_mutually_exclusive_group(required=True) # at leaset one of cif and xyz is provided
    # argparse will make sure only one of argument in structfile(xyz, cif) appear on command line
    structfile.add_argument("--xyz", type=str, default=None,
            help="The xyz structure file with second line specifying cell parameters")

    structfile.add_argument("--cif", type=str, default=None,
            help="The cif structure file")

    # run option
    subparser.add_argument("--runopt", type=str, default="gen",
            choices=["gen", "run", "genrun"],
            help="Generate or run or both at the same time.")

    subparser.add_argument("--auto", type=int, default=3,
            choices=[0, 1, 2, 3],
            help="auto:0 nothing, 1: copying files to server, 2: copying and executing, 3: pymatflow run inserver with direct submit,  in order use auto=1, 2, you must make sure there is a working ~/.pymatflow/server_[pbs|yh].conf")


    # -----------------------------------------------------------------
    #                      run params
    # -----------------------------------------------------------------

    subparser.add_argument("--mpi", type=str, default="",
            help="MPI command: like 'mpirun -np 4'")

    subparser.add_argument("--server", type=str, default="pbs",
            choices=["pbs", "yh"],
            help="type of remote server, can be pbs or yh")

    subparser.add_argument("--jobname", type=str, default="matflow-job",
            help="jobname on the pbs server")

    subparser.add_argument("--nodes", type=int, default=1,
            help="Nodes used in server")

    subparser.add_argument("--ppn", type=int, default=32,
            help="ppn of the server")




    # ==========================================================
    # transfer parameters from the arg subparser to static_run setting
    # ==========================================================

    args = parser.parse_args()

    # if no argument passed to matflow
    if len(sys.argv) == 1:
        # display help message when no args provided
        parser.print_help()
        sys.exit(1)

    # dealing wich structure files
    if args.xyz != None:
        xyzfile = args.xyz
    else:
        os.system("cif-to-xyz-modified.py -i %s -o %s.xyz" % (args.cif, args.cif))
        xyzfile = "%s.xyz" % args.cif




    if args.driver == "abinit":
        params = {}
        kpoints = {}

        params["chkprim"] = args.chkprim
        params["ecut"] = args.ecut
        params["ixc"] = args.ixc
        params["vdw_xc"] = args.vdw_xc
        params["vdw_tol"] = args.vdw_tol

        kpoints["kptopt"] = args.kptopt
        kpoints["ngkpt"] = args.ngkpt

        params["occopt"] = args.occopt
        params["nband"] = args.nband
        params["occ"] = args.occ

        params["nsppol"] = args.nsppol
        params["prtden"] = args.prtden
        params["prtdos"] = args.prtdos
        if args.runtype == 0:
            # static
            from pymatflow.abinit.static import static_run
            task = static_run()
            if get_kpath(args.kpath_manual, args.kpath_file) == None:
                print("================================================\n")
                print("Warning: matflow abinit\n")
                print("in abinit static runing you must provide kpath\n")
                sys.exit(1)
            task.dataset[3].electrons.kpoints.set_band(kptbounds=get_kpath(args.kpath_manual, args.kpath_file))
        elif args.runtype == 1:
            # optimization
            from pymatflow.abinit.opt import opt_run
            task = opt_run()
        elif args.runtype == 2:
            # dfpt-elastic-piezo-dielec
            from pymatflow.abinit.dfpt import dfpt_elastic_piezo_dielec
            task = dfpt_elastic_piezo_dielec()
        else:
            pass

        task.get_xyz(xyzfile)
        task.set_params(params=params)
        task.set_kpoints(kpoints=kpoints)
        task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
        task.run(directory=args.directory, runopt=args.runopt, auto=args.auto)


    elif args.driver == "cp2k":
        params = {}

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

        #force_eval["DFT-XC-VDW_POTENTIAL"] = args.usevdw
        params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-POTENTIAL_TYPE"] = args.vdw_potential_type
        params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-TYPE"] = args.pair_type
        params["FORCE_EVAL-DFT-XC-VDW_POTENTIAL-PAIR_POTENTIAL-R_CUTOFF"] = args.r_cutoff


        params["FORCE_EVAL-DFT-PRINT-ELF_CUBE-STRIDE"] = args.dft_print_elf_cube_stride
        params["FORCE_EVAL-DFT-PRINT-E_DENSITY_CUBE-STRIDE"] = args.dft_print_e_density_cube_stride


        params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-RANGE"] = args.properties_resp_slab_sampling_range
        params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-SURF_DIRECTION"] = args.properties_resp_slab_sampling_surf_direction
        params["FORCE_EVAL-PROPERTIES-RESP-SLAB_SAMPLING-ATOM_LIST"] = args.properties_resp_slab_sampling_atom_list
        if args.runtype == 0:
            from pymatflow.cp2k.static import static_run
            task = static_run()
            task.get_xyz(xyzfile)
            task.set_params(params=params)
            task.set_printout(option=args.printout_option)
            if 2 in args.printout_option and kpath != None:
                task.force_eval.dft.printout.band_structure.set_band(kpath=get_kpath(args.kpath_manual, args.kpath_file))
            task.set_vdw(usevdw=True if args.usevdw.lower() == "true" else False)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
            task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
        else:
            pass
    elif args.driver == "qe":
        control = {}
        electrons = {}
        system = {}
        ions = {}

        control["tstress"] = args.tstress
        system["ecutwfc"] = args.ecutwfc
        system["ecutrho"] = args.ecutrho
        system["occupations"] = args.occupations
        system["smearing"] = args.smearing
        system["degauss"] = args.degauss
        system["vdw_corr"] = args.vdw_corr
        system["nbnd"] = args.nbnd
        electrons["conv_thr"] = args.conv_thr

        system["nspin"] = args.nspin
        system["starting_magnetization"] = args.starting_magnetization
        system["noncolin"] = args.noncolin

        path = {}
        path["string_method"] = args.string_method
        path["nstep_path"] = args.nstep_path
        path["opt_scheme"] = args.opt_scheme
        path["num_of_images"] = args.num_of_images
        path["k_max"] = args.k_max
        path["k_min"] = args.k_min
        path["CI_scheme"] = args.ci_scheme
        path["path_thr"] = args.path_thr
        path["ds"] = args.ds
        path["first_last_opt"] = args.first_last_opt

        # for ph.x
        inputph = {}
        inputph["tr2_ph"] = args.tr2_ph
        inputph["lraman"] = args.lraman
        inputph["epsil"] = args.epsil
        inputph["nq1"] = args.nq[0]
        inputph["nq2"] = args.nq[1]
        inputph["nq3"] = args.nq[2]


        if args.runtype == 0:
            from pymatflow.qe.static import static_run
            task = static_run()
            task.get_xyz(xyzfile)
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_params(control=control, system=system, electrons=electrons)
            task.set_atomic_forces(pressure=args.pressure, pressuredir=args.pressuredir)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
            task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto)
        elif args.runtype == 1:
            pass
        elif args.runtype == 2:
            pass
        elif args.runtype == 3:
            from pymatflow.qe.neb import neb_run
            task = neb_run()
            task.get_images(images=args.images)
            task.set_kpoints(kpoints_option=args.kpoints_option, kpoints_mp=args.kpoints_mp)
            task.set_path(path=path)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
            task.neb(directory=directory, runopt=args.runopt, restart_mode=args.restart_mode, auto=args.auto)
        elif args.runtype == 4:
            from pymatflow.qe.dfpt import dfpt_run
            task = dfpt_run()
            task.get_xyz(xyzfile)
            task.set_inputph(inputph=inputph)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
            task.phx(directory=args.directory, runopt=args.runopt, auto=args.auto)

    elif args.driver == "siesta":
        params = {}

        params["MeshCutoff"] = args.meshcutoff
        params["SolutionMethod"] = args.solution_method
        params["XC.funtional"] = args.functional
        params["XC.authors"] = args.authors
        params["DM.Tolerance"] = args.tolerance
        params["DM.NumberPulay"] = args.numberpulay
        params["DM.MixingWeight"] = args.mixing
        params["OccupationFunction"] = args.occupation
        params["ElectronicTemperature"] = args.electronic_temperature

        if args.runtype == 0:
            from pymatflow.siesta.static import static_run
            task = static_run()
            task.get_xyz(xyzfile)

            task.properties.set_params(
                #bandlines = args.bandlines,
                #bandpoints = args.bandpoints,
                polarization_grids = args.polarization_grids,
                external_electric_field = args.external_electric_field,
                optical_energy_minimum = args.optical_energy_minimum,
                optical_energy_maximum = args.optical_energy_maximum,
                optical_broaden = args.optical_broaden,
                optical_scissor = args.optical_scissor,
                optical_mesh = args.optical_mesh,
                optical_polarization_type = args.optical_polarization_type,
                optical_vector = args.optical_vector,
                wannier90_unkgrid = args.wannier90_unkgrid,
                )

            if 3 in args.properties:
                task.properties.bandlines = bandlines

            task.set_params(params=params)
            task.set_kpoints(kpoints_mp=args.kpoints_mp)
            task.set_run(mpi=args.mpi, server=args.server, jobname=args.jobname, nodes=args.nodes, ppn=args.ppn)
            task.scf(directory=args.directory, runopt=args.runopt, auto=args.auto, properties=args.properties)

    # --------------------------------------------------------------------------



if __name__ == "__main__":
    main()
